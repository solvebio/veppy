import sys
from collections import namedtuple

from .features import EffectiveVariant
from .utils import codons
from .utils import fasta


Coordinate = \
    namedtuple(
        'Coordinate',
        [
            'start',
            'stop'
        ]
    )


CodingSequence = \
    namedtuple(
        'CodingSequence',
        [
            'coding_sequence',
            'genomic_coordinates',
            'coding_coordinates',
            'amino_acids',
            'amino_acid_numbers'
        ]
    )


BuilderResult = \
    namedtuple(
        'BuilderResult',
        [
            'coding_coordinates',
            'reference_sequence',
            'alternate_sequence',
            'upstream_sequence',
            'downstream_sequence'
        ]
    )


class CodingSequenceBuilder(object):
    @staticmethod
    def complement(x):
        return {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C',
        }.get(x, x)

    @staticmethod
    def complement_sequence(seq):
        # this will raise a KeyError if x isn't valid...
        return ''.join(
            reversed(
                [CodingSequenceBuilder.complement(x) for x in seq]
            )
        )

    def __init__(self, build, chromosome, transcript):
        self.fasta_file = fasta.GenbankFastaFile.for_build(build)
        self.transcript = transcript
        self.tree = self.transcript.coordinate_translator
        self.strand = self.transcript.strand

    def build_variant(self, orig_variant):

        # Helper functions for amino acids
        def normalize_coding_sequence(coding_coordinates, seq):
            """
            Remove "coding" nucleotides that are before/after coding region
            of transcript.
            """
            nucleotides = []
            coding_length = self.transcript.coding_offsets.coding_length
            for (coord, nucleotide) in zip(coding_coordinates, seq):
                if coord < 0 or coord > coding_length - 1:
                    nucleotides.append('.')
                else:
                    nucleotides.append(nucleotide)
            return ''.join(nucleotides)

        def get_three_letter_codes(coding_sequence):
            amino_acids = []
            for seq in codons.split_into_codons(coding_sequence):
                # even though we should **never** have a len(3) sequence in
                #  a Codon object that isn't an actual coding sequence, let's
                #  be safe :-)
                amino_acids.append(
                    (codons.get(seq) or {})
                    .get('three_letter_code', None))

            return amino_acids

        def get_amino_acids(coding_coordinates, coding_sequence):
            amino_acids = get_three_letter_codes(
                normalize_coding_sequence(
                    coding_coordinates, coding_sequence))

            # calculate amino acid numbers.
            # NOTE: 'coding_length' includes stop codon, so
            #  must subtract 3 to account for that b/c stop_codon
            #  has a 'null' amino acid number
            amino_acid_numbers = []
            coding_length = \
                self.transcript.coding_offsets.coding_length
            for i in range(len(coding_coordinates) / 3):
                coding_start = coding_coordinates[i * 3]
                if coding_start < 0 or \
                        (coding_start + 1) > (coding_length - 3):
                    # amino_acid_numbers.append(None)
                    amino_acid_numbers.append(-10000)
                else:
                    amino_acid_numbers.append(coding_start / 3 + 1)

            return (amino_acids, amino_acid_numbers)

        variant = EffectiveVariant(orig_variant)

        # need to clone original to get actual start/stop positions in
        #  event that start/stop get modified below when we're checking
        #  variant coordinates in relation to the transcript coding
        #  bounds. this is important in the case of DELETIONs that
        #  extend outside the transcript coding sequence
        variant_clone = variant.clone()

        # if transcript isn't coding or variant isn't in transcript
        #  coding region, return None
        if not self.transcript.is_coding or \
                variant.start > self.transcript.coding_stop or \
                variant.stop < self.transcript.coding_start:
            return None

        # bound variant! (work in genomic coordinates)
        if variant.stop > self.transcript.coding_stop:
            (lo, hi) = \
                (0, variant.stop - self.transcript.coding_stop)
            variant.stop = self.transcript.coding_stop
            variant.start = min(variant.start, variant.stop)
            variant.reference_allele = \
                variant.reference_allele[lo:hi]

        if variant.start < self.transcript.coding_start:
            (lo, hi) = \
                (self.transcript.coding_start - variant.start, sys.maxint)
            variant.start = self.transcript.coding_start
            variant.stop = max(variant.stop, variant.start)
            variant.reference_allele = \
                variant.reference_allele[lo:hi]

        alt_allele = variant.alternate_allele

        # -- Reference
        # build reference sequence and buffers
        #  symmetric codon buffer --> pull back N upstream and N downstream
        #   codons
        # NOTE: for insertions, [variant.start, variant.stop] overlaps the
        #  location and variant.stop is the location of the newly inserted
        #  allele(s)
        (
            coding_start,
            coding_stop,
            intron_start_offset,
            intron_stop_offset
        ) = self.tree.to_coding_range(variant.start, variant.stop)
        if variant.is_insertion:
            (coding_start, intron_start_offset) = \
                (coding_stop, intron_stop_offset)

        # TODO: how do we ever get here??
        if coding_start is None and coding_stop is None:
            return None

        # admittedly, a hack! if both intron start and stop offsets are None,
        #  variant is completely intronic
        # TODO: this is broken! what about exon deletion in which the deletion
        #  originates and terminates in separate introns????
        if intron_start_offset > 0 \
                and intron_stop_offset > 0:
            return None

        # "round" start and stop to nearest codon
        coding_start_round = coding_start / 3 * 3
        coding_stop_round = coding_stop / 3 * 3 + 2

        # HACK!
        # edge case: if variant is insertion (and indel??) and
        #  occurs as a codon boundary (between codons), then...
        (ref_coding_start, ref_coding_stop) = \
            (coding_start_round, coding_stop_round)
        if variant.is_insertion and (
            (self.strand == '+' and coding_start % 3 == 0) or
            (self.strand == '-' and coding_stop % 3 == 0)
        ):
            # NOTE: this '-1' is to account for the +1 offset when
            # calculating the downstream_coding ranges below
            coding_stop_round = coding_start_round - 1

            (ref_seq, ref_genomic_coordinates) = ('', [])
            ref_coding_coordinates = []
        else:

            (ref_seq, ref_genomic_coordinates) = self.build(
                variant.chromosome, ref_coding_start, ref_coding_stop)
            ref_coding_coordinates = \
                range(ref_coding_start, ref_coding_stop + 1)

        # build symmetric upstream/downstream codon buffer around ref
        # TODO: need to bound upstream_seq and downstream_seq
        _buffer = 3
        # upstream
        # if start > stop, build returns empty range
        (upstream_coding_start, upstream_coding_stop) = (
            coding_start_round - 3 * _buffer,
            coding_start_round - 1
        )

        (upstream_seq, upstream_genomic_coordinates) = \
            self.build(
                variant.chromosome, upstream_coding_start,
                upstream_coding_stop)
        upstream_coding_coordinates = \
            range(upstream_coding_start, upstream_coding_stop + 1)

        # downstream
        (downstream_coding_start, downstream_coding_stop) = (
            coding_stop_round + 1,
            coding_stop_round + 3 * _buffer
        )
        (downstream_seq, downstream_genomic_coordinates) = \
            self.build(
                variant.chromosome, downstream_coding_start,
                downstream_coding_stop)
        downstream_coding_coordinates = \
            range(downstream_coding_start, downstream_coding_stop + 1)

        # -- Alternate
        # build alt seq...
        # complement if negative strand...
        if self.strand == '-':
            (ref_seq, alt_allele, upstream_seq, downstream_seq) = map(
                CodingSequenceBuilder.complement_sequence,
                (ref_seq, alt_allele, upstream_seq, downstream_seq)
            )
            (
                ref_genomic_coordinates,
                upstream_genomic_coordinates,
                downstream_genomic_coordinates
            ) = map(
                lambda x: x[::-1],
                (
                    ref_genomic_coordinates,
                    upstream_genomic_coordinates,
                    downstream_genomic_coordinates
                )
            )

        if variant.is_insertion:
            i0 = coding_start - coding_start_round

            # if self.strand == '-':
            #     # edge case for insertions because of how they're encoded
            #     #  in Variant format...
            #     # TODO: better documentation
            #     i0 = coding_start - coding_start_round + 1

            i1 = i0

        else:
            i0 = coding_start - coding_start_round
            i1 = i0 + coding_stop - coding_start + 1

        alt_seq = ''.join([ref_seq[:i0], alt_allele, ref_seq[i1:]])

        # pad alt sequence to ensure it's has an integer number of codons
        _pad = len(alt_seq) % 3
        if _pad != 0:
            if self.strand == '+':
                alt_seq += self.fasta_file.get(
                    variant.chromosome,
                    variant_clone.stop + 1,
                    variant_clone.stop + (3 - _pad)
                )
            else:
                alt_seq += CodingSequenceBuilder.complement_sequence(
                    self.fasta_file.get(
                        variant.chromosome,
                        variant_clone.stop + 1,
                        variant_clone.stop + (3 - _pad)
                    )
                )

        # TODO: should we update ref_seq, alt_seq, etc. to remove
        #  non-coding nucleotides?
        (ref_amino_acids, ref_amino_acid_numbers) = \
            get_amino_acids(
                ref_coding_coordinates,
                ref_seq
        )

        alt_amino_acids = get_three_letter_codes(alt_seq)

        (upstream_amino_acids, upstream_amino_acid_numbers) = \
            get_amino_acids(
                upstream_coding_coordinates,
                upstream_seq
        )

        (downstream_amino_acids, downstream_amino_acid_numbers) = \
            get_amino_acids(
                downstream_coding_coordinates,
                downstream_seq
        )

        # NOTE: these are all CODING (!!!) sequences
        return BuilderResult(
            range(coding_start, coding_stop + 1),
            CodingSequence(
                ref_seq,
                ref_genomic_coordinates,
                ref_coding_coordinates,
                ref_amino_acids,
                ref_amino_acid_numbers
            ),
            CodingSequence(
                alt_seq,
                None,
                None,
                alt_amino_acids,
                None
            ),
            CodingSequence(
                upstream_seq,
                upstream_genomic_coordinates,
                upstream_coding_coordinates,
                upstream_amino_acids,
                upstream_amino_acid_numbers
            ),
            CodingSequence(
                downstream_seq,
                downstream_genomic_coordinates,
                downstream_coding_coordinates,
                downstream_amino_acids,
                downstream_amino_acid_numbers
            )
        )

    def build(self, chromosome, coding_start, coding_stop):
        # calculate genomic ranges
        genomic_ranges = self.tree.to_genomic_ranges(coding_start, coding_stop)

        # HACK: hmmmm.....
        if coding_start > coding_stop:
            return ('', [])

        # get genomic sequence
        # NOTE: genomic ranges are **ALWAYS** ordered from 5' to 3' w.r.t.
        #  the genomic coordinate system
        _seqs = []
        for _range in genomic_ranges:
            _seqs.append(
                self.fasta_file.get(
                    chromosome,
                    _range.start,
                    _range.stop
                )
            )
        seq = ''.join(_seqs)

        # output genomic coordinates of each bp in codon IN same order
        #  as codon
        codon_genomic_coordinates = []
        for gr in genomic_ranges:
            codon_genomic_coordinates.extend(
                range(gr.start, gr.stop + 1)
            )

        return (seq, codon_genomic_coordinates)
