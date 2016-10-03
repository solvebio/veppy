import sys
import re
import logging

from intervaltree import IntervalTree
from bisect import bisect
from collections import defaultdict, namedtuple

from . import errors

logger = logging.getLogger('veppy')


def intersects(start0, stop0, start1, stop1):
    return \
        overlaps(start0, stop0, start1, stop1) or \
        contains(start1, stop1, start0, stop0)


def overlaps(start0, stop0, start1, stop1):
    return start0 <= stop1 and stop0 >= start1


def contains(start0, stop0, start1, stop1):
    return start0 >= start1 and stop0 <= stop1


def contains_point(coord, start1, stop1):
    return contains(coord, coord, start1, stop1)


class Feature(object):
    def __init__(self, chromosome, start, stop, strand, frame,
                 info={}, genome_build='GRCh37'):
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.strand = strand
        self.frame = frame
        self.info = info
        self.genome_build = genome_build

    @property
    def length(self):
        return self.stop - self.start + 1

    def overlaps(self, feature):
        return overlaps(feature.start, feature.stop, self.start, self.stop)

    def contains(self, feature):
        return contains(feature.start, feature.stop, self.start, self.stop)

    def intersects(self, feature):
        return intersects(feature.start, feature.stop, self.start, self.stop)

    # override these methods to define custom sorting/comparator
    def __eq__(self, feature):
        return \
            self.__class__ == feature.__class__ and \
            self.chromosome == feature.chromosome and \
            self.start == feature.start and \
            self.stop == feature.stop

    def __ne__(self, obj):
        return not self.__eq__(obj)

    def __gt__(self, obj):
        return self.start > obj.start

    def __lt__(self, obj):
        return self.start < obj.start

    def __str__(self):
        return \
            'feature: %s, chromosome: %s, coordinates: [%s, %s]' % \
            (self.__class__.__name__, self.chromosome, self.start, self.stop)

    def __repr__(self):
        return str(self)

    def to_dict(self):
        return {
            'chromosome': self.chromosome,
            'start': self.start,
            'stop': self.stop,
            'strand': self.strand
        }


class VepSpliceAcceptorJunction(Feature):
    # defined in accordance with standard VEP definition:
    #  spans from [-8 bp, 3 bp] ([-3bp, 8bp]) of intron/exon (exon/intron)
    #  boundary
    def __init__(self, transcript, intron):
        chromosome = transcript.chromosome
        strand = transcript.strand
        frame = None

        if strand == '+':
            start = intron.stop - 7
            stop = intron.stop + 3

            self.splice_acceptor = \
                Feature(
                    chromosome,
                    intron.stop - 1,
                    intron.stop,
                    strand,
                    frame
                )
        else:
            start = intron.start - 3
            stop = intron.start + 7

            self.splice_acceptor = \
                Feature(
                    chromosome,
                    intron.start,
                    intron.start + 1,
                    strand,
                    frame
                )

        super(VepSpliceAcceptorJunction, self).__init__(
            chromosome, start, stop, strand, frame)


class VepSpliceDonorJunction(Feature):
    # defined in accordance with standard VEP definition:
    #  spans from [-8 bp, 3 bp] ([-3bp, 8bp]) of intron/exon (exon/intron)
    #  boundary
    def __init__(self, transcript, intron):
        chromosome = transcript.chromosome
        strand = transcript.strand
        frame = None

        if strand == '+':
            start = intron.start - 3
            stop = intron.start + 7

            self.splice_donor = \
                Feature(
                    chromosome,
                    intron.start,
                    intron.start + 1,
                    strand,
                    frame
                )
        else:
            start = intron.stop - 7
            stop = intron.stop + 3

            self.splice_donor = \
                Feature(
                    chromosome,
                    intron.stop - 1,
                    intron.stop,
                    strand,
                    frame
                )

        super(VepSpliceDonorJunction, self).__init__(
            chromosome, start, stop, strand, frame)


class ParentFeature(Feature):
    def __init__(self, *args, **kwargs):
        self.feature_trees = defaultdict(IntervalTree)
        self._built = False
        super(ParentFeature, self).__init__(*args, **kwargs)

    def to_dict(self):
        d = super(ParentFeature, self).to_dict()
        d['id'] = self.id
        return d

    @property
    def id(self):
        return self.info.get(self.id_field)

    # check ID for equality as well
    def __eq__(self, feature):
        return super(ParentFeature, self) == feature and \
            self.id == feature.id

    def has(self, feature_type):
        return len(self.feature_trees[feature_type]) > 0

    def first(self, feature_type):
        tree = self.feature_trees[feature_type]
        leaf = (list(tree[tree.begin()]) or [None])[0]
        return leaf.data if leaf else None

    def last(self, feature_type):
        tree = self.feature_trees[feature_type]
        leaf = (list(tree[tree.end() - 1]) or [None])[0]
        return leaf.data if leaf else None

    def add_feature(self, feature_type, feature):
        self.feature_trees[feature_type]\
            .addi(feature.start, feature.stop + 1, feature)

    def get_features(self, feature_type):
        # TODO: calling sorted every time we request a feature list
        #       could be pretty expensive...
        return sorted(
            [x.data for x in self.feature_trees[feature_type].items()]
        )

    def find_overlapping(self, feature_type, variant, first=True):
        start = variant.effective_variant.start
        stop = variant.effective_variant.stop
        features = \
            [l.data for l in
                self.feature_trees[feature_type.lower()][start:stop + 1]]

        if first:
            return (features or [None])[0]

        else:
            return features

    def find_overlapping2(self, feature_type, variant, first=True):
        # TODO: do we **really** just want to return the first?
        feature_array = self.feature_trees[feature_type.lower()]
        if len(feature_array) == 0:
            return None

        i = bisect(feature_array, variant)

        _features = []
        for j in {max(i - 1, 0), min(i, len(feature_array) - 1)}:
            feature = feature_array[j]
            if overlaps(
                variant.effective_variant.start,
                variant.effective_variant.stop,
                feature.start + (1 if variant.is_insertion else 0),
                feature.stop
            ):
                _features.append(feature)

            elif len(_features) > 0:
                break

        if first:
            return None if not _features else _features[0]

        else:
            return _features

    # wrapper method to safely build transcripts.
    #  e.g.  does not build correctly for SeqGene v.104 b/c of invalid,
    #   overlapping exons
    def build(self, data, force=False):
        if self._built and not force:
            return

        try:
            self._build(data, force=force)

            # assign numbers to exons, introns and cdses
            for feature_type in ['exon', 'intron']:
                features = self.get_features(feature_type)[::1 if self.strand == '+' else -1]  # noqa
                for (i, feature) in enumerate(features):
                    # start counting at 1!
                    feature.number = i + 1
                    feature.total = len(features)
        except Exception, e:
            logger.warn('Invalid feature %s' % self)
            logger.warn(e)
            self.errors = True

    def __str__(self):
        return \
            'feature: %s, %s, chromosome: %s, coordinates: [%s, %s]' % \
            (
                self.__class__.__name__,
                self.id,
                self.chromosome,
                self.start,
                self.stop
            )


class DefaultFeatureMixin(object):
    default = False


class Gene(DefaultFeatureMixin, ParentFeature):
    transcripts = None

    @property
    def symbol(self):
        return self.info.get(self.symbol_field)


class EnsemblGene(Gene):
    feature = 'gene'
    id_field = 'gene_id'
    symbol_field = 'gene_name'


class GencodeGene(EnsemblGene):
    pass


class NcbiGene(Gene):
    feature = 'GENE'
    id_field = 'gene_id'
    symbol_field = 'gene_name'


class Transcript(DefaultFeatureMixin, ParentFeature):
    errors = False
    protein_id = None

    CodingOffset = \
        namedtuple(
            'CodingOffset',
            [
                'coding_offset',
                'coding_length'
            ]
        )

    CodingRange = \
        namedtuple(
            'CodingRange',
            [
                'coding_start',
                'coding_stop',
                'intronic_offset_start',
                'intronic_offset_stop'
            ]
        )

    GenomicRange = \
        namedtuple(
            'GenomicRange',
            [
                'start',
                'stop'
            ]
        )

    @property
    def accession_number(self):
        # from ENST123456.78, return ENST123456
        return self.id.split('.', 1)[0]

    @property
    def accession_int(self):
        # from ENST123456.78, return 123456 as int
        # strips all non-digits from string and returns as int
        return int(re.sub(r"\D", "", self.id.split('.', 1)[0]))

    @property
    def accession_version(self):
        # from ENST123456.78, return 78
        if '.' in self.id:
            return int(self.id.rsplit('.', 1)[-1])
        return 0

    @property
    def gene_id(self):
        return self.info.get(self.gene_id_field)

    @property
    def exons(self):
        return self.get_features('exon')

    @property
    def coding_sequences(self):
        return self.get_features('cds')

    @property
    def untranslated_regions(self):
        return self.get_features('utr')

    @property
    def introns(self):
        return self.get_features('intron')

    @property
    def is_coding(self):
        return bool(self.coding_sequences)

    @property
    def coding_start(self):
        return self.first('cds').start

    @property
    def coding_stop(self):
        # TODO: fix this to use 'stop_codon' if it exists...
        return self.last('cds').stop

    @property
    def coding_offsets(self):
        return getattr(self, '_coding_offsets', None)

    # Coordinate Translator
    @property
    def coordinate_translator(self):
        if not hasattr(self, '_coordinate_translator'):
            self._coordinate_translator = \
                Transcript.CoordinateTranslator(
                    self.exons,
                    self.introns,
                    self.strand,
                    self.coding_offsets.coding_offset,
                    self.coding_offsets.coding_length
                )
        return self._coordinate_translator

    class CoordinateTranslator(object):
        class Leaf(object):
            def __init__(self, feature, coding_start, coding_stop):
                self.feature = feature
                self.start = feature.start
                self.stop = feature.stop
                self.coding_start = coding_start
                self.coding_stop = coding_stop

            def __str__(self):
                return 'genomic: [%s, %s], coding: [%s, %s]' % (
                    self.start, self.stop, self.coding_start, self.coding_stop
                )

        def __init__(
                self, exons, introns, strand, coding_offset, coding_length):
            self.strand = strand
            self.coding_offset = coding_offset
            self.coding_length = coding_length
            self._exon_tree = IntervalTree()
            self._intron_tree = IntervalTree()
            self._genomic_tree = IntervalTree()

            _coding_start = -self.coding_offset

            for exon in (exons if self.strand == '+' else exons[::-1]):
                leaf = Transcript.CoordinateTranslator.Leaf(
                    exon,
                    _coding_start,
                    _coding_start + exon.length - 1
                )

                self._genomic_tree.addi(leaf.start, leaf.stop + 1, leaf)
                self._exon_tree.addi(
                    leaf.coding_start, leaf.coding_stop + 1, leaf)

                # increment
                _coding_start = leaf.coding_stop + 1

            for intron in introns:
                # introns don't have coding coordinates, so use those of
                # adjacent exons
                leaf_genomic_upstream = \
                    list(self._genomic_tree[intron.start - 1])[0].data
                leaf_genomic_downstream = \
                    list(self._genomic_tree[intron.stop + 1])[0].data

                # NOTE: always assemble intronic offsets w.r.t. to the
                #  'coding stop' position of the upstream CDS
                if self.strand == '+':
                    leaf = \
                        Transcript.CoordinateTranslator.Leaf(
                            intron,
                            leaf_genomic_upstream.coding_stop,
                            leaf_genomic_downstream.coding_start
                        )
                else:
                    leaf = \
                        Transcript.CoordinateTranslator.Leaf(
                            intron,
                            leaf_genomic_downstream.coding_stop,
                            leaf_genomic_upstream.coding_start
                        )
                self._intron_tree.addi(leaf.start, leaf.stop + 1, leaf)

            # add introns that are upstream and downstream to the exon
            #  sequence
            # TODO: we may not need this, depending on how we choose to handle
            #  [start, stop] ranges that occur outside exon ranges
            if self.strand == '+':
                # straw upstream (genomic) intron
                straw0 = \
                    Feature('.', 0, self._genomic_tree.begin() - 1, self.strand, None)      # noqa
                leaf0 = \
                    Transcript.CoordinateTranslator.Leaf(straw0, -1, 0)
                self._intron_tree.addi(straw0.start, straw0.stop, leaf0)

                # straw downstream (genomic) intron
                straw1 = \
                    Feature('.', self._genomic_tree.end() + 1, sys.maxint, self.strand, None)      # noqa
                leaf1 = \
                    Transcript.CoordinateTranslator.Leaf(
                    straw1, self.coding_length - 1, self.coding_length)    # noqa
                self._intron_tree.addi(straw1.start, straw1.stop, leaf1)

            else:
                # straw upstream (genomic) intron
                straw0 = \
                    Feature('.', 0, self._genomic_tree.begin() - 1, self.strand, None)      # noqa
                leaf0 = \
                    Transcript.CoordinateTranslator.Leaf(straw0, self.coding_length - 1, self.coding_length)    # noqa

                self._intron_tree.addi(straw0.start, straw0.stop, leaf0)

                # straw downstream (genomic) intron
                straw1 = \
                    Feature('.', self._genomic_tree.end() + 1, sys.maxint, self.strand, None)      # noqa
                leaf1 = \
                    Transcript.CoordinateTranslator.Leaf(straw1, -1, 0)    # noqa
                self._intron_tree.addi(straw1.start, straw1.stop, leaf1)

        def to_coding_range(self, start, stop, hgvs_format=False):
            #  from above, introns have a coding_length == 1
            # TODO: set 'intron' attribute on leaves in '_intron_tree'
            #  above
            def _is_intron(leaf):
                return leaf.coding_stop - leaf.coding_start == 1

            # coding start
            range_coding_start = (
                list(self._genomic_tree[start] | self._intron_tree[start]) or
                [None]
            )[0]

            coding_start = None
            intron_coding_offset_start = 0
            leaf = range_coding_start.data
            if _is_intron(leaf):
                if self.strand == '+':
                    delta0 = start - leaf.start + 1
                    delta1 = leaf.stop + 1 - start
                    if hgvs_format and delta0 > delta1:
                        coding_start = leaf.coding_stop
                        intron_coding_offset_start = -delta1
                    else:
                        coding_start = leaf.coding_start
                        intron_coding_offset_start = delta0

                else:
                    delta0 = leaf.stop + 1 - stop
                    delta1 = stop - leaf.start + 1
                    if hgvs_format and delta0 > delta1:
                        coding_start = leaf.coding_stop
                        intron_coding_offset_start = -delta1
                    else:
                        coding_start = leaf.coding_start
                        intron_coding_offset_start = delta0
            else:
                if self.strand == '+':
                    coding_start = \
                        leaf.coding_start + (start - leaf.start)
                else:
                    coding_start = \
                        leaf.coding_start + (leaf.stop - stop)

            # coding stop
            range_coding_stop = (
                list(self._genomic_tree[stop] | self._intron_tree[stop]) or
                [None]
            )[0]

            coding_stop = None
            intron_coding_offset_stop = 0
            leaf = range_coding_stop.data
            if _is_intron(leaf):
                if self.strand == '+':
                    delta0 = stop - leaf.start + 1
                    delta1 = leaf.stop + 1 - stop
                    if hgvs_format and delta0 > delta1:
                        coding_stop = leaf.coding_stop
                        intron_coding_offset_stop = -delta1
                    else:
                        coding_stop = leaf.coding_start
                        intron_coding_offset_stop = delta0

                else:
                    delta0 = leaf.stop + 1 - start
                    delta1 = start - leaf.start + 1
                    if hgvs_format and delta0 > delta1:
                        coding_stop = leaf.coding_stop
                        intron_coding_offset_stop = -delta1
                    else:
                        coding_stop = leaf.coding_start
                        intron_coding_offset_stop = delta0

            else:
                if self.strand == '+':
                    coding_stop = \
                        leaf.coding_stop - (leaf.stop - stop)
                else:
                    coding_stop = \
                        leaf.coding_stop - (start - leaf.start)

            return \
                Transcript.CodingRange(
                    coding_start,
                    coding_stop,
                    intron_coding_offset_start,
                    intron_coding_offset_stop
                )

        def to_genomic_ranges(self, coding_start, coding_stop):
            genomic_ranges = []
            list_ranges = sorted(
                self._exon_tree[coding_start:coding_stop + 1],
                reverse=self.strand == '-'
            )

            for leaf in [r.data for r in list_ranges]:
                if self.strand == '+':
                    genomic_ranges.append(
                        Transcript.GenomicRange(
                            leaf.start + max(coding_start - leaf.coding_start, 0),  # noqa
                            leaf.stop - max(leaf.coding_stop - coding_stop, 0)   # noqa
                        )
                    )
                else:
                    genomic_ranges.append(
                        Transcript.GenomicRange(
                            leaf.start + max(leaf.coding_stop - coding_stop, 0),     # noqa
                            leaf.stop - max(coding_start - leaf.coding_start, 0)   # noqa
                        )
                    )

            return genomic_ranges

        def __str__(self):
            return 'coding sequences: %s' % map(str, self._tree)


class EnsemblTranscript(Transcript):
    feature = 'transcript'
    id_field = 'transcript_id'
    gene_id_field = 'gene_id'

    def _build(self, data, force=False):
        # create feature arrays...
        for d in data:
            feature_types = [d.get('feature').lower()]

            for feature_type in feature_types:
                self.add_feature(
                    feature_type,
                    Feature(
                        d['chromosome'],
                        d['start'],
                        d['stop'],
                        d['strand'],
                        d['frame'],
                        d['info']
                    )
                )

        # super hack! append 'stop_codon', if it exists, to
        #  the last CDS (if '+' strand) or first CDS (if '-' strand)
        #  including the stop codon in the last CDS makes lots of math
        #  and calculations easier for both VEP and HGVS
        stop_codon = (self.get_features('stop_codon') or [None])[-1]
        if self.strand == '+':
            stop_cds = self.last('cds')
        else:
            stop_cds = self.first('cds')

        if stop_codon and stop_cds:
            # BUG FIX: ENST00000470843
            #  transcript where entirety of last exon is a stop codon.
            #  if so, append a new CDS, else update stop position
            #  of last codon
            #  uggh....
            if self.strand == '+':
                self.add_feature(
                    'cds',
                    Feature(
                        stop_codon.chromosome,
                        stop_codon.start,
                        stop_codon.stop,
                        stop_codon.strand,
                        stop_codon.frame,
                        stop_codon.info
                    )
                )

            else:
                self.add_feature(
                    'cds',
                    Feature(
                        stop_codon.chromosome,
                        stop_codon.start,
                        stop_codon.stop,
                        stop_codon.strand,
                        stop_codon.frame,
                        stop_codon.info
                    )
                )

        # fix frame start...
        if self.strand == '+':
            start_cds = (self.coding_sequences or [None])[0]
            if start_cds and start_cds.frame is not None:
                start_cds.start += start_cds.frame
        else:
            start_cds = (self.coding_sequences or [None])[-1]
            if start_cds and start_cds.frame is not None:
                start_cds.stop -= start_cds.frame

        # calculate coding coordinate system
        if self.is_coding:
            if self.strand == '+':
                _cds = self.first('cds')
                _coding_offset = 0
                for _exon in self.exons:
                    if _exon.contains(_cds):
                        self._coding_offsets = \
                            Transcript.CodingOffset(
                                _coding_offset + (_cds.start - _exon.start),
                                sum([x.length for x in self.coding_sequences])
                            )
                        break

                    else:
                        _coding_offset += _exon.length

            else:
                _cds = self.last('cds')
                _coding_offset = 0
                for _exon in self.exons[::-1]:
                    if _exon.contains(_cds):
                        self._coding_offsets = \
                            Transcript.CodingOffset(
                                _coding_offset + (_exon.stop - _cds.stop),
                                sum([x.length for x in self.coding_sequences])
                            )
                        break

                    else:
                        _coding_offset += _exon.length

        # create introns...
        # not sure this sort is required anymore...
        sorted_exons = sorted(self.get_features('exon'))
        for i in range(1, len(sorted_exons)):
            x_prev = sorted_exons[i - 1]
            x = sorted_exons[i]
            intron = Feature(
                x.chromosome,
                x_prev.stop + 1,
                x.start - 1,
                x.strand,
                x.frame,
                x.info
            )

            self.add_feature('intron', intron)

            # build acceptor/donor splice junctions
            if self.strand == '+':
                self.add_feature(
                    'splice_acceptor_junction',
                    VepSpliceAcceptorJunction(self, intron)
                )

                self.add_feature(
                    'splice_donor_junction',
                    VepSpliceDonorJunction(self, intron)
                )

            else:
                # negative strand
                self.add_feature(
                    'splice_acceptor_junction',
                    VepSpliceAcceptorJunction(self, intron)
                )

                self.add_feature(
                    'splice_donor_junction',
                    VepSpliceDonorJunction(self, intron)
                )

        self._built = True


class GencodeTranscript(EnsemblTranscript):
    pass


class NcbiTranscript(Transcript):
    feature = 'RNA'
    id_field = 'transcript_id'
    gene_id_field = 'gene_id'

    def _build(self, data, force=False):
        # create feature arrays...
        for d in data:
            self.add_feature(
                d.get('feature').lower(),
                Feature(
                    d['chromosome'],
                    d['start'],
                    d['stop'],
                    d['strand'],
                    None,
                    d.get('info', {})
                )
            )

        # add in frames...
        # we have to calculate these by hand from the MapView file
        prev_cds = None
        for (i, cds) in enumerate(self.get_features('cds')):
            if i == 0:
                cds.frame = 0

            else:
                cds.frame = (prev_cds.frame + prev_cds.length) % 3

            prev_cds = cds

        # build exons
        exonic_features = \
            sorted(self.coding_sequences + self.get_features('utr'))
        exon = None
        for feature in exonic_features:
            if exon and exon.stop == feature.start - 1:
                exon.stop = feature.stop
                continue

            exon = \
                Feature(
                    feature.chromosome,
                    feature.start,
                    feature.stop,
                    feature.strand,
                    feature.frame,
                    feature.info
                )
            self.add_feature('exon', exon)

        # calculate coding coordinate system
        if self.is_coding:
            if self.strand == '+':
                _cds = (self.coding_sequences or [None])[0]
                _coding_offset = 0
                for _exon in self.exons:
                    if _exon.contains(_cds):
                        self._coding_offsets = \
                            Transcript.CodingOffset(
                                _coding_offset + (_cds.start - _exon.start),
                                sum([x.length for x in self.coding_sequences])
                            )
                        break

                    else:
                        _coding_offset += _exon.length

            else:
                _cds = (self.coding_sequences or [None])[-1]
                _coding_offset = 0
                for _exon in reversed(self.exons):
                    if _exon.contains(_cds):
                        self._coding_offsets = \
                            Transcript.CodingOffset(
                                _coding_offset + (_exon.stop - _cds.stop),
                                sum([x.length for x in self.coding_sequences])
                            )
                        break

                    else:
                        _coding_offset += _exon.length

        # create introns...
        # not sure this sort is required anymore...
        sorted_exons = sorted(self.get_features('exon'))
        for i in range(1, len(sorted_exons)):
            x_prev = sorted_exons[i - 1]
            x = sorted_exons[i]
            intron = Feature(
                x.chromosome,
                x_prev.stop + 1,
                x.start - 1,
                x.strand,
                x.frame,
                x.info
            )
            self.add_feature('intron', intron)

            # build acceptor/donor splice junctions
            if self.strand == '+':
                self.add_feature(
                    'splice_acceptor_junction',
                    VepSpliceAcceptorJunction(self, intron)
                )

                self.add_feature(
                    'splice_donor_junction',
                    VepSpliceDonorJunction(self, intron)
                )

            else:
                # negative strand
                self.add_feature(
                    'splice_acceptor_junction',
                    VepSpliceAcceptorJunction(self, intron)
                )

                self.add_feature(
                    'splice_donor_junction',
                    VepSpliceDonorJunction(self, intron)
                )
        self._built = True


class Variant(Feature):
    type = 'variant'

    def __init__(self, chromosome, start, ref, alt,
                 genome_build='GRCh37'):
        self.reference_allele = ref
        self.alternate_allele = alt

        super(Variant, self).__init__(
            chromosome,
            start,
            start + max(len(ref) - 1, 0),
            '.',
            None,
            genome_build=genome_build
        )
        self._offset_variant = None
        self._effective_variant = None

    def clone(self):
        return Variant(
            self.chromosome,
            self.start,
            self.reference_allele,
            self.alternate_allele,
            genome_build=self.genome_build
        )

    # override these methods to define custom sorting/comparator
    def hash(self):
        return str(self)

    @property
    def sequence_length(self):
        if self.is_insertion:
            return len(self.alternate_allele)

        else:
            return \
                abs(len(self.reference_allele) - len(self.alternate_allele))

    # Variant types
    @property
    def is_substitution(self):
        return \
            len(self.reference_allele) == len(self.alternate_allele)

    @property
    def is_indel(self):
        return \
            len(self.reference_allele) > 1 and len(self.alternate_allele) > 1

    @property
    def is_insertion(self):
        return \
            len(self.alternate_allele) > 1 and \
            len(self.reference_allele) == 1

    @property
    def is_deletion(self):
        return \
            len(self.reference_allele) > 1 and \
            len(self.alternate_allele) == 1

    @property
    def is_snp(self):
        return \
            len(self.reference_allele) == 1 and \
            len(self.reference_allele) == len(self.alternate_allele)

    @property
    def reference_length(self):
        return len(self.reference_allele)

    @property
    def alternate_length(self):
        return len(self.alternate_allele)

    # Variant sequence coordinates
    def offset_variant(self):
        if not self._offset_variant:
            variant = self.clone()

            if variant.reference_allele[0] == variant.alternate_allele[0]:
                variant.reference_allele = variant.reference_allele[1:]
                variant.alternate_allele = variant.alternate_allele[1:]
                variant.start += 1

            self._offset_variant = variant

        return self._offset_variant

    def offset_reference(self):
        if self.reference_allele[0] == self.alternate_allele[0]:
            return self.reference_allele[1:]

    def offset_alternate(self):
        if self.reference_allele[0] == self.alternate_allele[0]:
            return self.alternate_allele[1:]

    @property
    def effective_variant(self):
        if not self._effective_variant:
            self._effective_variant = EffectiveVariant(self)
        return self._effective_variant

    @property
    def sbid(self):
        return '{}-{}-{}-{}-{}'.format(
            self.genome_build,
            self.chromosome,
            self.start,
            self.stop,
            self.alternate_allele
        )

    def __eq__(self, feature):
        return \
            self.__class__ == feature.__class__ and \
            self.chromosome == feature.chromosome and \
            self.start == feature.start and \
            self.stop == feature.stop and \
            self.reference_allele == feature.reference_allele and \
            self.alternate_allele == feature.alternate_allele

    def __str__(self):
        return (
            'feature: %s, chromosome: %s, coordinates: [%s, %s, %s --> %s]' %     # noqa
            (
                self.__class__.__name__,
                self.chromosome,
                self.start,
                self.stop,
                self.reference_allele,
                self.alternate_allele
            )
        )


class EffectiveVariant(Variant):
    def __init__(self, variant):
        # avoid InceptionVariant. ensure that we don't make an
        #  EffectiveVariant for an EffectiveVariant
        if isinstance(variant, EffectiveVariant):
            chromosome = variant.chromosome
            start = variant.start
            reference_allele = variant.reference_allele
            alternate_allele = variant.alternate_allele

            # keep a reference to the original. useful for printing.
            self.original_variant = variant.original_variant
        else:
            # keep a reference to the original. useful for printing.
            self.original_variant = variant

            if variant.reference_allele == variant.alternate_allele:
                chromosome = variant.chromosome
                start = variant.start
                reference_allele = variant.reference_allele
                alternate_allele = variant.alternate_allele

            elif variant.is_snp:
                chromosome = variant.chromosome
                start = variant.start
                reference_allele = variant.reference_allele
                alternate_allele = variant.alternate_allele

            elif variant.is_indel:
                # remove common nucleotides...
                ref = variant.reference_allele
                alt = variant.alternate_allele
                i = 0
                while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
                    i += 1

                chromosome = variant.chromosome
                start = variant.start + i
                reference_allele = ref[i:]
                alternate_allele = alt[i:]

            elif variant.is_deletion:
                chromosome = variant.chromosome
                start = variant.start + 1
                reference_allele = variant.reference_allele[1:]
                alternate_allele = ''

            elif variant.is_insertion:
                chromosome = variant.chromosome
                start = variant.start
                reference_allele = ''
                alternate_allele = variant.alternate_allele[1:]

            else:
                raise errors.VeppyFeatureException(
                    'Unsupported variant: %s' % variant)

        super(EffectiveVariant, self).__init__(
            chromosome,
            start,
            reference_allele,
            alternate_allele
        )

        # force insertions to have a [start, stop] range that
        #  overlaps the location of the inserted alternate allele
        #  essentially, make variant.stop == site of first inserted
        #  alternate allele
        if self.is_insertion:
            self.stop = self.start + 1

    # Variant types
    @property
    def is_substitution(self):
        return self.original_variant.is_substitution

    @property
    def is_indel(self):
        return self.original_variant.is_indel

    @property
    def is_insertion(self):
        return self.original_variant.is_insertion

    @property
    def is_deletion(self):
        return self.original_variant.is_deletion

    @property
    def is_snp(self):
        return self.original_variant.is_snp

    def clone(self):
        return EffectiveVariant(self)

    def __str__(self):
        return str(self.original_variant)
