from .features import overlaps
from .features import contains
from .features import contains_point

from .utils.codons import synonymous
from .utils.codons import nonsynonymous
from .utils.codons import initiator_codon
from .utils.codons import stop_codon


def upstream_transcript_variant(variant, transcript):
    upstream_buffer = 5000
    if transcript.strand == '+':
        return \
            contains(
                variant.stop,
                variant.stop,
                transcript.start - 1 - (upstream_buffer - 1),
                transcript.start - 1,
            )
    else:
        return \
            contains(
                variant.start,
                variant.start,
                transcript.stop + 1,
                transcript.stop + 1 + (upstream_buffer - 1)
            )


def downstream_transcript_variant(variant, transcript):
    downstream_buffer = 2000
    if transcript.strand == '+':
        return \
            contains(
                variant.start,
                variant.start,
                transcript.stop + 1,
                transcript.stop + 1 + (downstream_buffer - 1)
            )
    else:
        return \
            contains(
                variant.stop,
                variant.stop,
                transcript.start - 1 - (downstream_buffer - 1),
                transcript.start - 1,
            )


# --------------------------------------------------
#  Splice Region, Donor and Acceptor
# --------------------------------------------------
def exon_loss(variant, transcript, exon):
    return variant.contains(exon)


# --------------------------------------------------
#  Splice Region, Donor and Acceptor
# --------------------------------------------------
def splice_donor_variant(variant, transcript, splice_model, *args):
    effective_variant = variant.effective_variant

    v_start = effective_variant.start
    v_stop = effective_variant.stop
    sd_start = splice_model.splice_donor.start
    sd_stop = splice_model.splice_donor.stop

    if variant.is_insertion:
        return v_stop == sd_stop
    else:
        return overlaps(v_start, v_stop, sd_start, sd_stop)


def splice_acceptor_variant(variant, transcript, splice_model, *args):
    effective_variant = variant.effective_variant

    v_start = effective_variant.start
    v_stop = effective_variant.stop
    sa_start = splice_model.splice_acceptor.start
    sa_stop = splice_model.splice_acceptor.stop

    if variant.is_insertion:
        return v_stop == sa_stop
    else:
        return overlaps(v_start, v_stop, sa_start, sa_stop)


def splice_region_variant(variant, transcript, splice_model, splice_region_type):  # noqa
    # genomic upstream
    if splice_region_type == 'acceptor':
        return _splice_acceptor_region_variant(
            variant, transcript, splice_model)

    else:
        return _splice_donor_region_variant(
            variant, transcript, splice_model)


def _splice_acceptor_region_variant(variant, transcript, splice_model):
    effective_variant = variant.effective_variant

    v_start = effective_variant.start
    v_stop = effective_variant.stop
    e_sr = splice_model.exon_splice_region
    i_sr = splice_model.intron_splice_region

    if transcript.strand == '+':
        if variant.is_insertion:
            return (
                contains_point(v_stop, i_sr.start + 1, i_sr.stop + 1) or
                contains_point(v_stop, e_sr.start, e_sr.stop)
            )
        else:
            return (
                overlaps(v_start, v_stop, i_sr.start, i_sr.stop) or
                overlaps(v_start, v_stop, e_sr.start, e_sr.stop)
            )

    else:
        if variant.is_insertion:
            return (
                contains_point(v_stop, i_sr.start, i_sr.stop) or
                contains_point(v_stop, e_sr.start + 1, e_sr.stop + 1)
            )
        else:
            return (
                overlaps(v_start, v_stop, i_sr.start, i_sr.stop) or
                overlaps(v_start, v_stop, e_sr.start, e_sr.stop)
            )


def _splice_donor_region_variant(variant, transcript, splice_model):
    effective_variant = variant.effective_variant

    v_start = effective_variant.start + \
        (1 if effective_variant.is_insertion else 0)
    v_stop = effective_variant.stop

    e_sr = splice_model.exon_splice_region
    i_sr = splice_model.intron_splice_region

    if transcript.strand == '+':
        if variant.is_insertion:
            return (
                contains_point(v_stop, i_sr.start, i_sr.stop) or
                contains_point(v_stop, e_sr.start + 1, e_sr.stop + 1)
            )
        else:
            return (
                overlaps(v_start, v_stop, i_sr.start, i_sr.stop) or
                overlaps(v_start, v_stop, e_sr.start, e_sr.stop)
            )

    else:
        if variant.is_insertion:
            return (
                contains_point(v_stop, i_sr.start + 1, i_sr.stop + 1) or
                contains_point(v_stop, e_sr.start, e_sr.stop)
            )
        else:
            return (
                overlaps(v_start, v_stop, i_sr.start, i_sr.stop) or
                overlaps(v_start, v_stop, e_sr.start, e_sr.stop)
            )

# --------------------------------------------------


def intronic_variant(variant, transcript, *args):
    return True


def five_prime_utr_variant(variant, transcript, utr):
    if not transcript.is_coding:
        return False

    if transcript.strand == '+':
        return utr.stop < transcript.coding_start

    else:
        return utr.start > transcript.coding_stop


def three_prime_utr_variant(variant, transcript, utr):
    if not transcript.is_coding:
        return False

    # BUG FIX: the +/- 3 was added to account for the fact that
    #  in some transcript libraries (GENCODE GTF) UTRs include the
    #  stop_codon whereas in others (ENCODE GTF) UTRs exclude the
    #  stop_codon
    if transcript.strand == '+':
        return (utr.start + 3) > transcript.coding_stop

    else:
        return (utr.stop - 3) < transcript.coding_start


def synonymous_variant(variant, transcript, ref_codon, variant_sequence):
    # BUG FIX: chr1, 201328373, G -> A, ENST00000438742
    #  this is a variant in the last base of an incompletely-annotated
    #  transcript. consequently, the len(ref_codon) < 3...
    if len(ref_codon) < 3:
        return False

    if not variant.is_snp and not variant.sequence_length % 3 == 0:
        return False

    return \
        ref_codon != variant_sequence and \
        variant.is_substitution and \
        synonymous(ref_codon, variant_sequence)


def missense_variant(variant, transcript, ref_codon, variant_sequence):
    # BUG FIX: chr1, 201328373, G -> A, ENST00000438742
    #  this is a variant in the last base of an incompletely-annotated
    #  transcript. consequently, the len(ref_codon) < 3...
    if len(ref_codon) < 3:
        return False

    if not variant.is_snp and not variant.sequence_length % 3 == 0:
        return False

    # check if a start...
    if initiator_codon_variant(
            variant, transcript, ref_codon, variant_sequence):
        return False

    return \
        variant.is_substitution and \
        nonsynonymous(ref_codon, variant_sequence) and \
        not stop_codon(ref_codon) and \
        not stop_codon(variant_sequence)


def frameshift_variant(variant, transcript, variant_sequence):
    return variant.sequence_length % 3 != 0


def inframe_insertion_variant(variant, transcript):
    return variant.is_insertion and variant.sequence_length % 3 == 0


def inframe_deletion_variant(variant, transcript):
    return variant.is_deletion and variant.sequence_length % 3 == 0


def inframe_indel_variant(variant, transcript):
    return (
        variant.is_indel and variant.sequence_length % 3 == 0 and
        len(variant.reference_allele) != len(variant.alternate_allele)
    )


# was ref mutated into a stop codon?
def stop_gained_variant(variant, transcript, ref_codon, variant_codon):
    return not stop_codon(ref_codon) and stop_codon(variant_codon)


# was ref mutated from a stop codon?
def stop_lost_variant(variant, transcript, ref_codon, variant_codon):
    return stop_codon(ref_codon) and not stop_codon(variant_codon)


# is synonymous mutation of stop codon?
def stop_retained_variant(variant, transcript, ref_codon, variant_codon):
    return \
        stop_codon(ref_codon) and synonymous(ref_codon, variant_codon)


# any change to first codon of a transcript...also has to be an ATG?
def initiator_codon_variant(variant, transcript, ref_codon, variant_codon):
    if transcript.strand == '+':
        return \
            initiator_codon(ref_codon) and \
            overlaps(
                variant.start,
                variant.stop,
                transcript.coding_start,
                transcript.coding_start + 2
            )

    else:
        return \
            initiator_codon(ref_codon) and \
            overlaps(
                variant.start,
                variant.stop,
                transcript.coding_stop - 2,
                transcript.coding_stop
            )


# was ref mutated from an initiator codon?
def start_lost_variant(variant, transcript, ref_codon, variant_codon):
    return initiator_codon_variant(
        variant, transcript, ref_codon, variant_codon
    ) and not initiator_codon(variant_codon)


# is synonymous mutation of an initiator codon?
def start_retained_variant(variant, transcript, ref_codon, variant_codon):
    return initiator_codon_variant(
        variant, transcript, ref_codon, variant_codon
    ) and initiator_codon(variant_codon)


# transcript with an incomplete terminal codon...
def incomplete_terminal_codon_variant(
    variant,
    transcript,
    ref_codon,
    variant_codon
):
    if transcript.strand == '+':
        return \
            ref_codon != variant_codon and \
            overlaps(
                variant.start,
                variant.stop,
                transcript.coding_stop - 2,
                transcript.coding_stop
            ) and \
            not stop_codon(ref_codon)

    else:
        return \
            ref_codon != variant_codon and \
            overlaps(
                variant.start,
                variant.stop,
                transcript.coding_start,
                transcript.coding_start + 2
            ) and \
            not stop_codon(ref_codon)


def feature_elongation_variant(variant, transcript):
    return variant.is_insertion or (
        variant.is_indel and
        len(variant.reference_allele) < len(variant.alternate_allele)
    )


def feature_truncation_variant(variant, transcript):
    return variant.is_deletion or (
        variant.is_indel and
        len(variant.reference_allele) > len(variant.alternate_allele)
    )


def nmd_transcript_variant(variant, transcript, cds, **kwargs):
    if transcript.strand == '+':
        if cds == transcript.coding_sequences[-1]:
            return True
        elif len(transcript.coding_sequences) > 1:
            return overlaps(
                variant.start,
                variant.stop,
                transcript.coding_sequences[-2].stop - 44,
                transcript.coding_sequences[-2].stop
            )
        else:
            return False

    else:
        if cds == transcript.coding_sequences[0]:
            return True
        elif len(transcript.coding_sequences) > 1:
            return overlaps(
                variant.start,
                variant.stop,
                transcript.coding_sequences[1].start + 44,
                transcript.coding_sequences[1].start
            )
        else:
            return False


def non_coding_transcript_variant(variant, transcript, *args, **kwargs):
    return \
        transcript.overlaps(variant) and \
        not transcript.has('cds')


# redundant call to "find_overlapping('exon'...)"
def non_coding_transcript_exon_variant(variant, transcript):
    return \
        not transcript.has('cds') and \
        transcript.find_overlapping('exon', variant) is not None


def coding_sequence_variant(variant, transcript, *args):
    return True


def no_sequence_alteration(variant, *args):
    return variant.reference_allele == variant.alternate_allele


# def mature_miRNA_variant(variant, transcript, distance=0):
#     raise NotImplementedError()


# def regulatory_region_variant(variant, transcript, distance=0):
#     raise NotImplementedError()


# def TF_binding_site_variant(variant, transcript, distance=0):
#     raise NotImplementedError()


# def transcript_ablation(variant, transcript, distance=0):
#     raise NotImplementedError()


# def transcript_amplification(variant, transcript, distance=0):
#     raise NotImplementedError()


# def TFBS_ablation(variant, transcript, distance=0):
#     raise NotImplementedError()


# def TFBS_amplification(variant, transcript, distance=0):
#     raise NotImplementedError()


# def regulatory_region_ablation(variant, transcript, distance=0):
#     raise NotImplementedError()


# def regulatory_region_amplification(variant, transcript, distance=0):
#     raise NotImplementedError()


# def transcript_elongation(variant, transcript, distance=0):
#     raise NotImplementedError()


# def transcript_truncation(variant, transcript, distance=0):
#     raise NotImplementedError()
