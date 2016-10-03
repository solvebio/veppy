import os
import logging
from glob import glob
from itertools import islice
from abc import ABCMeta
from time import time
from collections import defaultdict
from intervaltree import IntervalTree

from .features import EnsemblGene
from .features import EnsemblTranscript
from .features import GencodeGene
from .features import GencodeTranscript
from .features import NcbiGene
from .features import NcbiTranscript

from .utils import readers
from .utils import helpers

from . import errors
from . import SUPPORTED_BUILDS
from . import FEATURE_LIMIT


logger = logging.getLogger('veppy')


class FeatureFile(object):
    __metaclass__ = ABCMeta

    class BoundedCache(object):
        def __init__(self, size=1024):
            self._cache = [None] * size

        def key(self, feature):
            return hash(feature.id) % len(self._cache)

        def contains(self, feature):
            _feature = self.get(feature)
            if _feature:
                return _feature.id == feature.id
            return False

        def get(self, feature):
            return self._cache[self.key(feature)]

        def put(self, feature):
            self._cache[self.key(feature)] = feature

    def __init__(self, filepath, reindex=False):
        self._filepath = os.path.abspath(filepath)
        self._inspect_indexes()
        self.index(reindex=reindex)
        self._gene_trees = None
        self._genes_map = None
        self._transcript_trees = None
        self._transcripts_map = None
        self._protein_accessions_map = None

        self._cache = FeatureFile.BoundedCache()

    def genes(self, chromosome):
        self._load_transcripts()
        return self._gene_trees[chromosome]

    def transcripts(self, chromosome):
        self._load_transcripts()
        return self._transcript_trees[chromosome]

    def get_gene(self, chromosome, gene_id, build=True):
        self._load_transcripts()
        return self._genes_map.get(chromosome).get(gene_id)

    def get_transcript(self, transcript_id, build=True):
        self._load_transcripts()

        t = self._transcripts_map.get(transcript_id)
        if t and build:
            self.build_feature(t)

        # only return Transcripts that do not have build errors
        #  (see 'build_feature' below)
        if t and not t.errors:
            return t

        else:
            return None

    def get_transcript_by_protein_accession(self, protein_accession,
                                            build=True):
        self._load_transcripts()

        t = self._protein_accessions_map.get(protein_accession)
        if t and build:
            self.build_feature(t)

        # only return Transcripts that do not have build errors
        #  (see 'build_feature' below)
        if t and not t.errors:
            return t

        else:
            return None

    # get all overlapping features...
    def find(
            self,
            features,
            variant,
            walk_distance=None,
            upstream_buffer=0,
            downstream_buffer=0,
            build=True):
        # ensure effective variant!
        chromosome = variant.effective_variant.chromosome
        start = variant.effective_variant.start
        stop = variant.effective_variant.stop

        # add upstream/downstream buffers
        # ensure valid interval.
        #  if adjusted start > adjusted stop, swap
        start += upstream_buffer
        stop -= downstream_buffer
        if stop < start:
            (start, stop) = (stop, start)

        return \
            self.find_overlapping(
                features, chromosome, start, stop, build=build)

    def find_overlapping(self, features, chromosome, start, stop, build=True):
        _features = [r.data for r in features(chromosome)[start:stop + 1]]

        if build:
            for f in _features:
                if self._cache.contains(f):
                    f = self._cache.get(f)
                else:
                    self.build_feature(f)
                    self._cache.put(f)

        return [_f for _f in _features if not _f.errors]

    def find_genes(self, variant, upstream_buffer=0, downstream_buffer=0,
                   build=True):
        return self.find(
            self.genes,
            variant,
            upstream_buffer=upstream_buffer,
            downstream_buffer=downstream_buffer,
            build=build
        )

    def find_transcripts(self, variant, upstream_buffer=0, downstream_buffer=0,
                         build=True):
        transcripts = self.find(
            self.transcripts,
            variant,
            upstream_buffer=upstream_buffer,
            downstream_buffer=downstream_buffer,
            build=build
        )

        return transcripts

    def load(self, feature_id):
        _fp = os.path.join(
            self.parent_feature_index_directory,
            feature_id + '.' + self.extension
        )

        if not os.path.exists(_fp):
            logger.warn(
                'feature \'%s\' not found in index.' % feature_id)
            return []

        with self.reader_class(_fp) as r:
            return [data for (line, data) in r]

    def index(self, reindex=False):
        if not reindex and self._indexed:
            return

        # create index...
        # on-disk index structure:
        #  <root folder>/<transcript_index_filepath>
        #  <root folder>/<parent_feature>/<parent_feature_id>.<extension>
        # where 'root folder' is the basename of the source filepath minus
        #  the extension
        def _flush(parent_id, batch):
            if not parent_id or not batch:
                return

            _filepath = os.path.join(
                self.parent_feature_index_directory,
                parent_id + '.' + self.extension
            )
            if not os.path.exists(os.path.dirname(_filepath)):
                os.makedirs(os.path.dirname(_filepath))

            with open(_filepath, 'w') as f_out:
                f_out.write('\n'.join(batch) + '\n')

        _genes = []
        _transcripts = []
        _proteins = []
        _transcript_id = None
        _transcript_subfeatures = None
        _transcript_protein = None

        logger.info('Indexing \'%s\'' % self._filepath)
        t0 = time()
        for (line, data) in \
                islice(self.reader_class(self._filepath), FEATURE_LIMIT):
            if not data:
                continue

            if data['feature'] == self.transcript_class.feature:
                _flush(_transcript_id, _transcript_subfeatures)
                _transcript_subfeatures = []
                _transcript_id = data['info'][self.transcript_class.id_field]
                _transcripts.append(line)

                # protein to transcript mapping lines
                if _transcript_protein:
                    _proteins.append(_transcript_protein)
                _transcript_protein = None

                if len(_transcripts) % 1000 == 0:
                    logger.info(
                        '...indexed %s transcripts' % len(_transcripts)
                    )

            elif data['feature'] == self.gene_class.feature:
                # for the time being, we're not including any "subfeatures"
                #  in a gene model
                _genes.append(line)
                if len(_genes) % 1000 == 0:
                    logger.info('...indexed %s genes' % len(_genes))

            else:
                if _transcript_id:
                    if data['feature'] == 'CDS' and \
                            _transcript_protein is None:
                        _transcript_protein = line
                    _transcript_subfeatures.append(line)

        # flush last transcript batch...
        _flush(_transcript_id, _transcript_subfeatures)

        if _genes:
            with open(self.gene_index_file, 'w') as f_out:
                f_out.write('\n'.join(_genes) + '\n')

        if _transcripts:
            with open(self.transcript_index_file, 'w') as f_out:
                f_out.write('\n'.join(_transcripts) + '\n')

        if _proteins:
            with open(self.protein_index_file, 'w') as f_out:
                f_out.write('\n'.join(_proteins) + '\n')

        self._indexed = True

        logger.info(
            'Done! Indexed %s genes, %s transcripts in %ss.' %
            (len(_genes), len(_transcripts), int(time() - t0))
        )

    def build_feature(self, feature):
        feature.build(self.load(feature.id))

    @classmethod
    def for_build(klass, build):
        translated_build = helpers.supported_build(build)
        if not translated_build:
            raise errors.VeppyFileException(
                'Unsupported build \'%s\'.' % build)
        if not klass.BUILD_CACHE.get(translated_build):
            klass.BUILD_CACHE[translated_build] = \
                klass(klass.filepath_for_build(translated_build))

        return klass.BUILD_CACHE[translated_build]

    @classmethod
    def filepath_for_build(klass, build):
        (release, version) = helpers.parse_build(build)
        logger.info(
            'Loading build %s (release: %s, version: %s).' %
            (build, release, version or 'None')
        )

        return helpers.feature_filepath(build, klass.gene_model)

    @property
    def gene_index_file(self):
        (p, x) = os.path.splitext(self._filepath)
        return os.path.join(
            p,
            os.path.basename(p) + '.%s.%s' %
            (self.gene_class.feature, x.lstrip('.'))
        )

    @property
    def transcript_index_file(self):
        (p, x) = os.path.splitext(self._filepath)
        return os.path.join(
            p,
            os.path.basename(p) + '.%s.%s' %
            (self.transcript_class.feature, x.lstrip('.'))
        )

    @property
    def protein_index_file(self):
        (p, x) = os.path.splitext(self._filepath)
        return os.path.join(
            p,
            os.path.basename(p) + '.%s.%s' %
            ('protein', x.lstrip('.'))
        )

    @property
    def parent_feature_index_directory(self):
        return os.path.join(
            os.path.splitext(self._filepath)[0],
            self.transcript_class.feature
        )

    def _inspect_indexes(self):
        if \
                os.path.exists(self.gene_index_file) and \
                os.path.exists(self.protein_index_file) and \
                os.path.exists(self.transcript_index_file) and \
                os.path.exists(self.parent_feature_index_directory) and \
                glob(
                    os.path.join(
                        self.parent_feature_index_directory.rstrip('/'),
                        '*.' + self.extension
                    )
                ):
            self._indexed = True

        else:
            self._indexed = False

    # ensure genes loaded first!
    #  transcript loading links gene model to transcript model
    def _load_genes(self):
        if not self._indexed:
            self.index()

        if all(
            [
                self._gene_trees,
                self._genes_map
            ]
        ):
            return

        self._genes = defaultdict(list)
        self._gene_trees = defaultdict(IntervalTree)
        self._genes_map = defaultdict(dict)
        for (line, data) in \
                islice(self.reader_class(self.gene_index_file), FEATURE_LIMIT):       # noqa

            chromosome = data['chromosome']
            gene = self.gene_class(
                chromosome,
                data['start'],
                data['stop'],
                data['strand'],
                data['frame'],
                data['info']
            )
            self._gene_trees[chromosome].addi(
                gene.start, gene.stop + 1, gene)
            self._genes_map[chromosome][gene.id] = gene

    def _load_transcripts(self):
        if not self._indexed:
            self.index()

        # ensure genes are loaded!
        self._load_genes()
        if all(
            [
                self._transcript_trees,
                self._transcripts_map,
            ]
        ):
            return

        self._transcripts = defaultdict(list)
        self._transcript_trees = defaultdict(IntervalTree)
        self._transcripts_map = dict()
        self._protein_accessions_map = dict()
        cnt = 0
        for (line, data) in \
                islice(self.reader_class(self.transcript_index_file), FEATURE_LIMIT):       # noqa

            chromosome = data['chromosome']
            transcript = self.transcript_class(
                chromosome,
                data['start'],
                data['stop'],
                data['strand'],
                data['frame'],
                data['info']
            )
            self._transcript_trees[chromosome].addi(
                transcript.start, transcript.stop + 1, transcript)
            self._transcripts_map[transcript.id] = transcript
            cnt += 1

            # get associated Gene by ID
            gene = self.get_gene(transcript.chromosome, transcript.gene_id)
            # attach Gene model to Transcript
            # and Transcript to Gene model
            transcript.gene = gene
            if gene.transcripts is None:
                gene.transcripts = []
            gene.transcripts.append(transcript)

        # now, load proteins...
        self._load_proteins()

    def _load_proteins(self):
        if not self._indexed:
            self.index()

        # ensure transcripts are loaded!
        self._load_transcripts()
        if all(
            [
                self._protein_accessions_map
            ]
        ):
            return

        for (line, data) in \
                islice(self.reader_class(self.protein_index_file), FEATURE_LIMIT):       # noqa

            protein_id = data['info']['protein_id']
            transcript_id = data['info']['transcript_id']
            transcript = self.get_transcript(transcript_id, build=False)
            transcript.protein_id = protein_id
            self._protein_accessions_map[protein_id] = transcript


# ENSEMBL GTF

class ExtendedGtfReader(readers.GtfReader):
    def _process_next(self, line):
        return (line, super(ExtendedGtfReader, self)._process_next(line))


class GtfFeatureFile(FeatureFile):
    extension = 'gtf'
    reader_class = ExtendedGtfReader


class EnsemblGtfFile(GtfFeatureFile):
    BUILD_CACHE = {k: None for k in SUPPORTED_BUILDS}

    gene_model = 'ensembl'
    gene_class = EnsemblGene
    transcript_class = EnsemblTranscript


# GENCODE GTF

class GencodeGtfReader(ExtendedGtfReader):
    def _process_next(self, line):
        (line, _next) = super(GencodeGtfReader, self)._process_next(line)
        # clean chromosome!
        if _next:
            _next['chromosome'] = _next['chromosome'].lstrip('cChHrR')
            if 'protein_id' not in _next['info']:
                _next['info']['protein_id'] = None

        return (line, _next)


class GencodeGtfFile(EnsemblGtfFile):
    BUILD_CACHE = {k: None for k in SUPPORTED_BUILDS}

    gene_model = 'gencode'
    gene_class = GencodeGene
    transcript_class = GencodeTranscript
    reader_class = GencodeGtfReader


# NCBI MapView

class ExtendedNcbiMapViewReader(readers.NcbiMapViewReader):
    def _process_next(self, line):
        _next = super(ExtendedNcbiMapViewReader, self)._process_next(line)
        # keep GENEs and 'NM_' transcripts...
        if _next['feature_type'] != 'GENE' and \
                not _next['transcript'].startswith('NM_'):
            return None

        return (
            line,
            {
                'feature': _next['feature_type'],
                'chromosome': _next['chromosome'],
                'start': _next['chr_start'],
                'stop': _next['chr_stop'],
                'strand': _next['chr_orient'],
                'frame': None,          # MapView file doesn't contain frame...
                'info': {
                    'transcript_id': _next['transcript'],
                    'gene_id': _next['feature_id'],
                    'gene_name': _next['feature_name'] if _next['feature_type'] == 'GENE' else None,       # noqa
                    'protein_id': _next['feature_name'] if _next['feature_type'] == 'CDS' else None       # noqa
                }
            }
        )


class NcbiMapViewFile(FeatureFile):
    BUILD_CACHE = {k: None for k in SUPPORTED_BUILDS}

    gene_model = 'ncbi'
    extension = 'md'
    reader_class = ExtendedNcbiMapViewReader
    gene_class = NcbiGene
    transcript_class = NcbiTranscript


MODEL_MAP = {
    'gencode': GencodeGtfFile,
    'ncbi': NcbiMapViewFile,
    'refseq': NcbiMapViewFile
}


def get_model(gene_model):
    if gene_model.lower() not in MODEL_MAP:
        raise errors.VeppyFileException(
            'Unsupported gene model \'%s\'.' % gene_model)

    return MODEL_MAP[gene_model.lower()]
