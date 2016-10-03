import functions


def to_dict(csq):
    return {
        'so_accession': csq['SO_accession'],
        'so_term': csq['SO_term'],
        'impact': csq['impact'],
        'rank': csq['rank']
    }


class Consequence(object):
    @staticmethod
    def impact_score(imp):
        return \
            {
                impacts.HIGH: 10 ** 6,
                impacts.MODERATE: 10 ** 5,
                impacts.LOW: 10 ** 4,
                impacts.MODIFIER: 10 ** 3,
                impacts.NONE: 10 ** 0
            }[imp]

    def __init__(self, consequence):
        self._dict = to_dict(consequence)
        self.so_accession = self._dict['so_accession']
        self.so_term = self._dict['so_term']
        self.impact = self._dict['impact']
        self.rank = self._dict['rank']
        self._score = \
            Consequence.impact_score(self.impact) - \
            self.rank

    def __hash__(self):
        return hash(self.so_term)

    def __eq__(self, obj):
        return self.so_term == obj.so_term

    def __ne__(self, obj):
        return not self.__eq__(obj)

    def __gt__(self, obj):
        if self._score == obj._score:
            return self.so_term > obj.so_term

        return self._score > obj._score

    def __lt__(self, obj):
        if self._score == obj._score:
            return self.so_term < obj.so_term

        return self._score < obj._score

    def __str__(self):
        return '%s, %s' % \
            (self.__class__.__name__, self._dict)


class impacts(object):
    HIGH = 'HIGH'
    MODERATE = 'MODERATE'
    LOW = 'LOW'
    MODIFIER = 'MODIFIER'
    NONE = 'NONE'


# impact source: http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
# --------------------------------------------------
#  impact: HIGH
# --------------------------------------------------
exon_loss = {
    'SO_accession': 'SO:0001572',
    'SO_term': 'exon_loss',
    'description': 'A sequence variant whereby an exon is lost from the transcript.',  # noqa
    'function': functions.exon_loss,
    'impact': impacts.HIGH,
    'rank': 2,
}

frameshift_variant = {
    'SO_accession': 'SO:0001589',
    'SO_term': 'frameshift_variant',
    'description': 'A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three',  # noqa
    'function': functions.frameshift_variant,
    'impact': impacts.HIGH,
    'rank': 3,
}

splice_acceptor_variant = {
    'SO_accession': 'SO:0001574',
    'SO_term': 'splice_acceptor_variant',
    'description': 'A splice variant that changes the 2 base region at the 3\' end of an intron',  # noqa
    'function': functions.splice_acceptor_variant,
    'impact': impacts.HIGH,
    'rank': 7,
}

splice_donor_variant = {
    'SO_accession': 'SO:0001575',
    'SO_term': 'splice_donor_variant',
    'description': 'A splice variant that changes the 2 base region at the 5\' end of an intron',  # noqa
    'function': functions.splice_donor_variant,
    'impact': impacts.HIGH,
    'rank': 8,
}

start_lost_variant = {
    'SO_accession': 'SO:0002012',
    'SO_term': 'start_lost',
    'description': 'A codon variant that changes at least one base of the canonical start codon.',  # noqa
    'function': functions.start_lost_variant,
    'impact': impacts.HIGH,
    'rank': 6,
}

stop_gained_variant = {
    'SO_accession': 'SO:0001587',
    'SO_term': 'stop_gained',
    'description': 'A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript',  # noqa
    'function': functions.stop_gained_variant,
    'impact': impacts.HIGH,
    'rank': 4,
}

stop_lost_variant = {
    'SO_accession': 'SO:0001578',
    'SO_term': 'stop_lost',
    'description': 'A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript',  # noqa
    'function': functions.stop_lost_variant,
    'impact': impacts.HIGH,
    'rank': 5,
}

# --------------------------------------------------
#  impact: MODERATE
# --------------------------------------------------
inframe_insertion_variant = {
    'SO_accession': 'SO:0001821',
    'SO_term': 'inframe_insertion',
    'description': 'An inframe non synonymous variant that inserts bases into in the coding sequence',  # noqa
    'function': functions.inframe_insertion_variant,
    'impact': impacts.MODERATE,
    'rank': 11,
}

inframe_deletion_variant = {
    'SO_accession': 'SO:0001822',
    'SO_term': 'inframe_deletion',
    'description': 'An inframe non synonymous variant that deletes bases from the coding sequence',  # noqa
    'function': functions.inframe_deletion_variant,
    'impact': impacts.MODERATE,
    'rank': 13,
}

inframe_indel_variant = {
    'SO_accession': 'SO:0001820',
    'SO_term': 'inframe_indel',
    'description': 'A coding sequence variant where the change does not alter the frame of the transcript.',  # noqa
    'function': functions.inframe_indel_variant,
    'impact': impacts.MODERATE,
    'rank': 12,
}

missense_variant = {
    'SO_accession': 'SO:0001583',
    'SO_term': 'missense_variant',
    'description': 'A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved',  # noqa
    'function': functions.missense_variant,
    'impact': impacts.MODERATE,
    'rank': 10,
}

splice_region_variant = {
    'SO_accession': 'SO:0001630',
    'SO_term': 'splice_region_variant',
    'description': 'A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron',  # noqa
    'function': functions.splice_region_variant,
    'impact': impacts.MODERATE,
    'rank': 18,
}

# --------------------------------------------------
#  impact: LOW
# --------------------------------------------------
incomplete_terminal_codon_variant = {
    'SO_accession': 'SO:0001626',
    'SO_term': 'incomplete_terminal_codon_variant',
    'description': 'A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed',  # noqa
    'function': functions.incomplete_terminal_codon_variant,
    'impact': impacts.LOW,
    'rank': 19,
}

initiator_codon_variant = {
    'SO_accession': 'SO:0001582',
    'SO_term': 'initiator_codon_variant',
    'description': 'A codon variant that changes at least one base of the first codon of a transcript',  # noqa
    'function': functions.initiator_codon_variant,
    'impact': impacts.LOW,
    'rank': 21,
}

start_retained_variant = {
    'SO_accession': 'SO:0002019',
    'SO_term': 'start_retained_variant',
    'description': 'A sequence variant where at least one base in the start codon is changed, but the start remains.',  # noqa
    'function': functions.start_retained_variant,
    'impact': impacts.LOW,
    'rank': 20,
}

stop_retained_variant = {
    'SO_accession': 'SO:0001567',
    'SO_term': 'stop_retained_variant',
    'description': 'A sequence variant where at least one base in the terminator codon is changed, but the terminator remains',  # noqa
    'function': functions.stop_retained_variant,
    'impact': impacts.LOW,
    'rank': 24,
}

synonymous_variant = {
    'SO_accession': 'SO:0001819',
    'SO_term': 'synonymous_variant',
    'description': 'A sequence variant where there is no resulting change to the encoded amino acid',  # noqa
    'function': functions.synonymous_variant,
    'impact': impacts.LOW,
    'rank': 22,
}

# --------------------------------------------------
#  impact: MODIFIER
# --------------------------------------------------
five_prime_utr_variant = {
    'SO_accession': 'SO:0001623',
    'SO_term': '5_prime_UTR_variant',
    'description': 'A UTR variant of the 5\' UTR',
    'function': functions.five_prime_utr_variant,
    'impact': impacts.MODIFIER,
    'rank': 26,
}

three_prime_utr_variant = {
    'SO_accession': 'SO:0001624',
    'SO_term': '3_prime_UTR_variant',
    'description': 'A UTR variant of the 3\' UTR',
    'function': functions.three_prime_utr_variant,
    'impact': impacts.MODIFIER,
    'rank': 27,
}

coding_sequence_variant = {
    'SO_accession': 'SO:0001580',
    'SO_term': 'coding_sequence_variant',
    'description': 'A sequence variant that changes the coding sequence',
    'function': functions.coding_sequence_variant,
    'impact': impacts.MODIFIER,     # NOTE: this is also listed as MODERATE in source   # noqa
    'rank': 41,     # NOTE: this also has a rank of 25
}

# downstream_gene_variant = {
#     'SO_accession': 'SO:0001632',
#     'SO_term': 'downstream_gene_variant',
#     'description': 'A sequence variant located 3\' of a gene',
#     'function': functions.downstream_gene_variant,
#     'impact': impacts.MODIFIER,
#     'rank': 30,
# }

downstream_transcript_variant = {
    'SO_accession': 'SO:0001987',
    'SO_term': 'downstream_transcript_variant',
    'description': 'A feature variant, where the alteration occurs downstream of the transcript TSS.',       # noqa
    'function': functions.downstream_transcript_variant,
    'impact': impacts.MODIFIER,
    'rank': 30,
}

feature_elongation_variant = {
    'SO_accession': 'SO:0001907',
    'SO_term': 'feature_elongation',
    'description': 'A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence',  # noqa
    'function': functions.feature_elongation_variant,
    'impact': impacts.MODIFIER,
    'rank': 39,
}

feature_truncation_variant = {
    'SO_accession': 'SO:0001906',
    'SO_term': 'feature_truncation',
    'description': 'A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence',  # noqa
    'function': functions.feature_truncation_variant,
    'impact': impacts.MODIFIER,
    'rank': 38,
}

intergenic_variant = {
    'SO_accession': 'SO:0001628',
    'SO_term': 'intergenic_variant',
    'description': 'A sequence variant located in the intergenic region, between genes',  # noqa
    'function': lambda x: True,
    'impact': impacts.MODIFIER,
    'rank': 40,
}

intronic_variant = {
    'SO_accession': 'SO:0001627',
    'SO_term': 'intron_variant',
    'description': 'A transcript variant occurring within an intron',
    'function': functions.intronic_variant,
    'impact': impacts.MODIFIER,
    'rank': 37,
}

nmd_transcript_variant = {
    'SO_accession': 'SO:0001621',
    'SO_term': 'nmd_transcript_variant',
    'description': 'A variant in a transcript that is the target of NMD',  # noqa
    'function': functions.nmd_transcript_variant,
    'impact': impacts.MODIFIER,
    'rank': 25,
}

non_coding_transcript_exon_variant = {
    'SO_accession': 'SO:0001792',
    'SO_term': 'non_coding_exon_variant',
    'description': 'A sequence variant that changes non-coding exon sequence in a non-coding transcript',  # noqa
    'function': functions.non_coding_transcript_exon_variant,
    'impact': impacts.MODIFIER,
    'rank': 42,
}

non_coding_transcript_variant = {
    'SO_accession': 'SO:0001619',
    'SO_term': 'non_coding_transcript_variant',
    'description': 'A transcript variant of a non coding RNA gene',
    'function': functions.non_coding_transcript_variant,
    'impact': impacts.MODIFIER,
    'rank': 43,
}

upstream_transcript_variant = {
    'SO_accession': 'SO:0001986',
    'SO_term': 'upstream_transcript_variant',
    'description': 'A feature variant, where the alteration occurs upstream of the transcript TSS.',       # noqa
    'function': functions.upstream_transcript_variant,
    'impact': impacts.MODIFIER,
    'rank': 29,
}

# upstream_gene_variant = {
#     'SO_accession': 'SO:0001631',
#     'SO_term': 'upstream_gene_variant',
#     'description': 'A sequence variant located 5\' of a gene',
#     'function': functions.upstream_gene_variant,
#     'impact': impacts.MODIFIER,
#     'rank': 29,
# }

# --------------------------------------------------
#  impact: NONE
# --------------------------------------------------
no_sequence_alteration = {
    'SO_accession': '',
    'SO_term': 'no_sequence_alteration',
    'description': 'A variant with identical reference and alternate alleles.',
    'function': functions.no_sequence_alteration,
    'impact': impacts.NONE,
    'rank': 100,
}

# Null variant
NULL_VARIANT_CONSEQUENCES = set([
    Consequence(stop_gained_variant),
    Consequence(frameshift_variant),
    Consequence(splice_acceptor_variant),
    Consequence(splice_donor_variant),
    Consequence(start_lost_variant),
    Consequence(exon_loss)
])


def is_null(csq):
    return csq in NULL_VARIANT_CONSEQUENCES


# Util for testing..
_term_map = {}
for csq in [
    frameshift_variant,
    splice_acceptor_variant,
    splice_donor_variant,
    start_lost_variant,
    stop_gained_variant,
    stop_lost_variant,
    inframe_insertion_variant,
    inframe_deletion_variant,
    missense_variant,
    splice_region_variant,
    incomplete_terminal_codon_variant,
    initiator_codon_variant,
    start_retained_variant,
    stop_retained_variant,
    synonymous_variant,
    five_prime_utr_variant,
    three_prime_utr_variant,
    coding_sequence_variant,
    downstream_transcript_variant,
    # downstream_gene_variant,
    feature_elongation_variant,
    feature_truncation_variant,
    # intergenic_variant,
    intronic_variant,
    non_coding_transcript_exon_variant,
    non_coding_transcript_variant,
    upstream_transcript_variant
    # upstream_gene_variant
]:
    _term_map[csq['SO_term']] = csq
