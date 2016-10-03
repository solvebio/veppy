import logging

from . import features
from . import errors
from . import consequences
from .feature_file import get_model
from .utils import helpers
from .sequence import CodingSequenceBuilder
from .splice_model import SpliceJunctionModel
from .utils.codons import split_into_codons

logger = logging.getLogger('veppy')


# wrapper to enable special-casing of insertions
def find_overlapping(transcript, feature_type, variant, first=True):
    overlapping_features = transcript.find_overlapping(
        feature_type, variant, first=first)

    # insertions are constructed such that [start, stop] contains
    #  the location of the inserted nucleotides, e.g.:
    #   (start, stop, allele) = (10, 11, 'ACT') means that 'ACT' was
    #   inserted immediately after position 10
    if not variant.is_insertion:
        return overlapping_features

    if feature_type == 'exon':
        return overlapping_features

    elif feature_type == 'cds':
        return overlapping_features

    elif feature_type == 'utr':
        return overlapping_features

    elif feature_type == 'intron':
        return overlapping_features

    elif feature_type == 'splice_acceptor_junction':
        return overlapping_features

    elif feature_type == 'splice_donor_junction':
        return overlapping_features

    raise features.VeppyFeatureException(
        'Unknown feature type \'%s\'' % feature_type)


class EffectPredictionResults(object):
    def __init__(self, variant, results=None):
        self.variant = variant

        if results:
            self.results = results
        else:
            self.results = []


class EffectPredictionModel(object):
    splice_model_class = SpliceJunctionModel

    def __init__(
        self,
        genome_build,
        transcript,
        variant
    ):
        self.genome_build = genome_build
        self.transcript = transcript
        self.variant = variant
        self.consequences = []
        self.context = {}

    def result(self, most_severe=False):
        _result = {
            'transcript': self.transcript,
            'consequences': self.consequences
        }
        # update results model with context fields...
        for k in [
            'cdses',
            'exons',
            'introns',
            'sequences',
            'coding_coordinates'
        ]:
            _result[k] = self.context.get(k, None)

        # sort consequences so that most severe effect is first
        if len(_result['consequences']) > 0:
            _sorted = sorted(_result['consequences'], reverse=True)
            if most_severe:
                _result['consequences'] = [_sorted[0]]

            else:
                _result['consequences'] = _sorted

        return _result

    # mutates result 'context'
    def update_context(self, **kwargs):
        self.context.update(kwargs)

    # mutates csqs
    def add_result(self, consequence, result, stop=False):
        if not result:
            return

        self.consequences.append(
            consequences.Consequence(consequence))

        if stop:
            raise errors.StopEffectPrediction()

    def calculate_splicing_effects(self):
        transcript = self.transcript
        variant = self.variant
        splice_donor_junction = \
            find_overlapping(transcript, 'splice_donor_junction', variant)

        if splice_donor_junction:
            splice_model = self.splice_model_class.get_donor_junction_model(
                transcript, splice_donor_junction)

            # check splice donor!
            for c in [
                consequences.splice_donor_variant,
                consequences.splice_region_variant
            ]:
                result = \
                    c.get('function')(
                        variant,
                        transcript,
                        splice_model,
                        'donor'
                    )
                self.add_result(c, result)

        splice_acceptor_junction = \
            find_overlapping(transcript, 'splice_acceptor_junction', variant)

        if splice_acceptor_junction:
            splice_model = self.splice_model_class.\
                get_acceptor_junction_model(
                    transcript, splice_acceptor_junction)

            # check splice acceptor!
            for c in [
                consequences.splice_acceptor_variant,
                consequences.splice_region_variant
            ]:
                result = \
                    c.get('function')(
                        variant,
                        transcript,
                        splice_model,
                        'acceptor'
                    )
                self.add_result(c, result)

    def calculate_nmd_effect(self):
        transcript = self.transcript
        variant = self.variant
        cds = self.context.get('cdses', [None])[0]
        if not cds:
            return

        csqs = self.consequences
        nmd_csq = consequences.nmd_transcript_variant
        for csq in csqs:
            if consequences.is_null(csq) and \
                    nmd_csq['function'](variant, transcript, cds):
                # we're done
                self.consequences.append(consequences.Consequence(nmd_csq))         # noqa
                break

    def calculate_consequences(self):
        build = self.genome_build
        transcript = self.transcript
        variant = self.variant

        # always check!
        # ref == alt?
        for c in [consequences.no_sequence_alteration]:
            result = c.get('function')(variant, transcript)
            self.add_result(c, result, stop=True)

        # non-coding transcripts don't have features,
        for c in [consequences.non_coding_transcript_variant,
                  consequences.non_coding_transcript_exon_variant]:
            result = c.get('function')(variant, transcript)
            self.add_result(c, result)

        # upstream/downstream?
        for c in [consequences.upstream_transcript_variant,
                  consequences.downstream_transcript_variant]:
            result = c.get('function')(variant, transcript)

            # we're done if is upstream or downstream
            self.add_result(c, result, stop=True)

        # exon loss trumps **ALL**
        exons = find_overlapping(transcript, 'exon', variant, first=False)
        if exons:
            self.update_context(exons=exons)
        for exon in exons:
            for c in [consequences.exon_loss]:
                result = c.get('function')(variant, transcript, exon)
                self.add_result(c, result, stop=True)

        # utr
        utr = find_overlapping(transcript, 'utr', variant)
        if utr:
            for c in [consequences.five_prime_utr_variant,
                      consequences.three_prime_utr_variant]:
                result = c.get('function')(variant, transcript, utr)
                self.add_result(c, result)

        # elongation/truncation
        for c in [consequences.feature_elongation_variant,
                  consequences.feature_truncation_variant]:
            result = c.get('function')(variant, transcript)
            self.add_result(c, result)

        # get cds
        cds = find_overlapping(transcript, 'cds', variant)
        if cds:
            self.update_context(cdses=[cds])

        # is intronic??
        intron = find_overlapping(transcript, 'intron', variant)
        if intron:
            self.update_context(introns=[intron])
            for c in [consequences.intronic_variant]:
                result = c.get('function')(variant, transcript)
                self.add_result(c, result)

        # Splice Region, Donor, and Acceptor!
        # splice_donor_or_acceptor = False
        self.calculate_splicing_effects()

        # non-coding, nothing left to do...
        if not transcript.is_coding:
            raise errors.StopEffectPrediction()

        # elif splice_donor_or_acceptor:
        #     # BUG FIX: insertions on the (genomic coordinate) upstream
        #     #  side of a CDS will not "overlap" an intron, even though
        #     #  we may call such a variant a "splice acceptor/donor"
        #     #  for consistency, if splice donor/acceptor, variant also
        #     #  must be classified as intronic
        #     self.add_result(consequences.intronic_variant, True)

        # check if this is a coding variant...
        if not cds:
            # not coding, so punt! we've already checked 3'/5' and
            #  upstream/downstream variants
            raise errors.StopEffectPrediction()

        if transcript.strand == '+':
            coding_variant = cds.overlaps(variant) and \
                not (variant.start >= transcript.coding_stop - 2)
        else:
            coding_variant = cds.overlaps(variant) and \
                not (variant.stop <= transcript.coding_start + 2)

        # coding sequence calculation and effect prediction
        # TODO: check that transcript is cached in feature file.
        builder = \
            CodingSequenceBuilder(build, variant.chromosome, transcript)
        builder_result = builder.build_variant(variant)

        if not builder_result:
            logger.debug('Intronic variant encountered %s.' % variant)
            return

        # build codons and save coding sequences in features model
        codons_ref_seq = split_into_codons(
            builder_result.reference_sequence.coding_sequence)
        codons_var_seq = split_into_codons(
            builder_result.alternate_sequence.coding_sequence)

        self.update_context(
            coding_coordinates=builder_result.coding_coordinates,   # noqa
            sequences={
                'reference_sequence': builder_result.reference_sequence,
                'alternate_sequence': builder_result.alternate_sequence,
                'upstream_sequence': builder_result.upstream_sequence,
                'downstream_sequence': builder_result.downstream_sequence
            }
        )

        # check terminal (initiator/stop) effects first
        # initiator/start
        for c in [consequences.initiator_codon_variant,
                  consequences.start_retained_variant,
                  consequences.start_lost_variant]:
            for codon_ref_seq in codons_ref_seq:
                result = c.get('function')(
                    variant,
                    transcript,
                    codon_ref_seq,
                    codons_var_seq[0]
                )
                self.add_result(
                    c, result, stop=(c == consequences.start_lost_variant))

        # stop (call stop_gained below)
        for c in [consequences.incomplete_terminal_codon_variant,
                  consequences.stop_retained_variant,
                  consequences.stop_lost_variant]:
            for codon_ref_seq in codons_ref_seq:
                result = c.get('function')(
                    variant,
                    transcript,
                    codon_ref_seq,
                    codons_var_seq[0]
                )
                self.add_result(
                    c, result, stop=(c == consequences.stop_lost_variant))

        # check frameshift/inframe insertions and deletions
        if coding_variant:
            for c in [consequences.frameshift_variant]:
                result = \
                    c.get('function')(
                        variant,
                        transcript,
                        None
                    )
                # frameshift punt. if frameshift, we're done!
                self.add_result(c, result, stop=True)

            for c in [consequences.inframe_deletion_variant,
                      consequences.inframe_insertion_variant,
                      consequences.inframe_indel_variant]:
                result = \
                    c.get('function')(
                        variant,
                        transcript
                    )
                self.add_result(c, result)

        # Coding functions
        # check first codon of var seq and ref seq for
        # NOTE: only checking the first codon has the important (intended)
        #  effect of possibly causing us to miss (or incorrectly handle)
        #  M --> N (M > 1, N > 1) INDELS. e.g.:
        # chromosome: 1, start: 35351702, GGG -> GTT
        # NOTE: added logic to support calling missense, synonymous variants
        #  for ALL codons in an N --> N variant
        if len(codons_ref_seq) == len(codons_var_seq):
            num_codons = len(codons_ref_seq)
        else:
            num_codons = 1

        # synonymous/missense
        for c in [consequences.coding_sequence_variant,
                  consequences.synonymous_variant,
                  consequences.missense_variant]:
            for i in range(0, num_codons):
                codon_ref = codons_ref_seq[i]
                codon_var = codons_var_seq[i]
                result = \
                    c.get('function')(
                        variant,
                        transcript,
                        codon_ref,
                        codon_var
                    )
                self.add_result(c, result)

        # **ALSO**, check var sequence for stop gained
        for c in [consequences.stop_gained_variant]:
            for codon_var_seq in codons_var_seq:
                result = \
                    c.get('function')(
                        variant,
                        transcript,
                        codons_ref_seq[0],
                        codon_var_seq
                    )
                self.add_result(c, result)


def calculate_consequences(genome_build, chromosome, start, reference, allele,
                           most_severe=False, default_only=False,
                           gene_model='ncbi',
                           effect_prediction_model=EffectPredictionModel,
                           exclude_neighbors=False):
    """
    The main entrypoint to veppy.
    Returns predicted effects for a given variant.
    """

    feature_file = get_model(gene_model).for_build(genome_build)

    # TODO: handle this case more reasonably.
    #  Possibly iterate through 'A', 'C', 'G', 'T' if alt == 'N'
    ref = reference.upper()
    alt = allele.upper()

    # construct effective Variant...
    variant_orig = features.Variant(chromosome, start, ref, alt)
    variant = features.EffectiveVariant(variant_orig)

    if exclude_neighbors:
        transcripts = feature_file.find_transcripts(
            variant, upstream_buffer=0, downstream_buffer=0)
    else:
        transcripts = feature_file.find_transcripts(
            variant, upstream_buffer=5000, downstream_buffer=2000)

    # intergenic
    if len(transcripts) == 0:
        return EffectPredictionResults(variant_orig)

    helpers.set_defaults(transcripts)
    if default_only:
        transcripts = [t for t in transcripts if t.default and t.gene.default]
        if not transcripts:
            logger.warn(
                'Attempted to calculate effects for default transcript '
                'on default gene only, but none found!')
            return EffectPredictionResults(variant_orig)

    # 'results' dict gets updated by 'consequences_for_model'
    ep_models = []
    for transcript in transcripts:
        ep_models.append(
            effect_prediction_model(
                genome_build,
                transcript,
                variant
            )
        )

        # short circuit if alt == 'N' or ref == alt
        if alt in {'N'}:
            continue

        try:
            ep_models[-1].calculate_consequences()
        except errors.StopEffectPrediction:
            pass

        # check nmd_transcript_variant
        #  NMD only applies to null variants (defined in consequences.py)
        ep_models[-1].calculate_nmd_effect()

    return EffectPredictionResults(
        variant_orig,
        [m.result(most_severe=most_severe) for m in ep_models]
    )
