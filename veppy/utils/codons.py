from collections import defaultdict

CODONS_FULL = {
    # Alanine
    'GCA': {
        'one_letter_code': 'A',
        'three_letter_code': 'Ala',
        'amino_acid': 'Alanine'
    },
    'GCC': {
        'one_letter_code': 'A',
        'three_letter_code': 'Ala',
        'amino_acid': 'Alanine'
    },
    'GCG': {
        'one_letter_code': 'A',
        'three_letter_code': 'Ala',
        'amino_acid': 'Alanine'
    },
    'GCT': {
        'one_letter_code': 'A',
        'three_letter_code': 'Ala',
        'amino_acid': 'Alanine'
    },

    # Cysteine
    'TGC': {
        'one_letter_code': 'C',
        'three_letter_code': 'Cys',
        'amino_acid': 'Cysteine'
    },
    'TGT': {
        'one_letter_code': 'C',
        'three_letter_code': 'Cys',
        'amino_acid': 'Cysteine'
    },

    # Aspartic Acid
    'GAC': {
        'one_letter_code': 'D',
        'three_letter_code': 'Asp',
        'amino_acid': 'Aspartic Acid'
    },
    'GAT': {
        'one_letter_code': 'D',
        'three_letter_code': 'Asp',
        'amino_acid': 'Aspartic Acid'
    },

    # Glutamic Acid
    'GAA': {
        'one_letter_code': 'E',
        'three_letter_code': 'Glu',
        'amino_acid': 'Glutamic acid'
    },
    'GAG': {
        'one_letter_code': 'E',
        'three_letter_code': 'Glu',
        'amino_acid': 'Glutamic acid'
    },

    # Phenylalanine
    'TTC': {
        'one_letter_code': 'F',
        'three_letter_code': 'Phe',
        'amino_acid': 'Phenylalanine'
    },
    'TTT': {
        'one_letter_code': 'F',
        'three_letter_code': 'Phe',
        'amino_acid': 'Phenylalanine'
    },

    # Glycine
    'GGA': {
        'one_letter_code': 'G',
        'three_letter_code': 'Gly',
        'amino_acid': 'Glycine'
    },
    'GGC': {
        'one_letter_code': 'G',
        'three_letter_code': 'Gly',
        'amino_acid': 'Glycine'
    },
    'GGG': {
        'one_letter_code': 'G',
        'three_letter_code': 'Gly',
        'amino_acid': 'Glycine'
    },
    'GGT': {
        'one_letter_code': 'G',
        'three_letter_code': 'Gly',
        'amino_acid': 'Glycine'
    },

    # Histidine
    'CAC': {
        'one_letter_code': 'H',
        'three_letter_code': 'His',
        'amino_acid': 'Histidine'
    },
    'CAT': {
        'one_letter_code': 'H',
        'three_letter_code': 'His',
        'amino_acid': 'Histidine'
    },

    # Isoleucine
    'ATA': {
        'one_letter_code': 'I',
        'three_letter_code': 'Ile',
        'amino_acid': 'Isoleucine'
    },
    'ATC': {
        'one_letter_code': 'I',
        'three_letter_code': 'Ile',
        'amino_acid': 'Isoleucine'
    },
    'ATT': {
        'one_letter_code': 'I',
        'three_letter_code': 'Ile',
        'amino_acid': 'Isoleucine'
    },

    # Lysine
    'AAA': {
        'one_letter_code': 'K',
        'three_letter_code': 'Lys',
        'amino_acid': 'Lysine'
    },
    'AAG': {
        'one_letter_code': 'K',
        'three_letter_code': 'Lys',
        'amino_acid': 'Lysine'
    },

    # Leucine
    'CTA': {
        'one_letter_code': 'L',
        'three_letter_code': 'Leu',
        'amino_acid': 'Leucine'
    },
    'CTC': {
        'one_letter_code': 'L',
        'three_letter_code': 'Leu',
        'amino_acid': 'Leucine'
    },
    'CTG': {
        'one_letter_code': 'L',
        'three_letter_code': 'Leu',
        'amino_acid': 'Leucine'
    },
    'CTT': {
        'one_letter_code': 'L',
        'three_letter_code': 'Leu',
        'amino_acid': 'Leucine'
    },
    'TTA': {
        'one_letter_code': 'L',
        'three_letter_code': 'Leu',
        'amino_acid': 'Leucine'
    },
    'TTG': {
        'one_letter_code': 'L',
        'three_letter_code': 'Leu',
        'amino_acid': 'Leucine'
    },

    # Methionine (Start Codon)
    'ATG': {
        'one_letter_code': 'M',
        'three_letter_code': 'Met',
        'amino_acid': 'Methionine'
    },

    # Asparagine
    'AAC': {
        'one_letter_code': 'N',
        'three_letter_code': 'Asn',
        'amino_acid': 'Asparagine'
    },
    'AAT': {
        'one_letter_code': 'N',
        'three_letter_code': 'Asn',
        'amino_acid': 'Asparagine'
    },

    # Proline
    'CCA': {
        'one_letter_code': 'P',
        'three_letter_code': 'Pro',
        'amino_acid': 'Proline'
    },
    'CCC': {
        'one_letter_code': 'P',
        'three_letter_code': 'Pro',
        'amino_acid': 'Proline'
    },
    'CCG': {
        'one_letter_code': 'P',
        'three_letter_code': 'Pro',
        'amino_acid': 'Proline'
    },
    'CCT': {
        'one_letter_code': 'P',
        'three_letter_code': 'Pro',
        'amino_acid': 'Proline'
    },

    # Glutamine
    'CAA': {
        'one_letter_code': 'Q',
        'three_letter_code': 'Gln',
        'amino_acid': 'Glutamine'
    },
    'CAG': {
        'one_letter_code': 'Q',
        'three_letter_code': 'Gln',
        'amino_acid': 'Glutamine'
    },

    # Arginine
    'AGA': {
        'one_letter_code': 'R',
        'three_letter_code': 'Arg',
        'amino_acid': 'Arginine'
    },
    'AGG': {
        'one_letter_code': 'R',
        'three_letter_code': 'Arg',
        'amino_acid': 'Arginine'
    },
    'CGA': {
        'one_letter_code': 'R',
        'three_letter_code': 'Arg',
        'amino_acid': 'Arginine'
    },
    'CGC': {
        'one_letter_code': 'R',
        'three_letter_code': 'Arg',
        'amino_acid': 'Arginine'
    },
    'CGG': {
        'one_letter_code': 'R',
        'three_letter_code': 'Arg',
        'amino_acid': 'Arginine'
    },
    'CGT': {
        'one_letter_code': 'R',
        'three_letter_code': 'Arg',
        'amino_acid': 'Arginine'
    },

    # Serine
    'AGC': {
        'one_letter_code': 'S',
        'three_letter_code': 'Ser',
        'amino_acid': 'Serine'
    },
    'AGT': {
        'one_letter_code': 'S',
        'three_letter_code': 'Ser',
        'amino_acid': 'Serine'
    },
    'TCA': {
        'one_letter_code': 'S',
        'three_letter_code': 'Ser',
        'amino_acid': 'Serine'
    },
    'TCC': {
        'one_letter_code': 'S',
        'three_letter_code': 'Ser',
        'amino_acid': 'Serine'
    },
    'TCG': {
        'one_letter_code': 'S',
        'three_letter_code': 'Ser',
        'amino_acid': 'Serine'
    },
    'TCT': {
        'one_letter_code': 'S',
        'three_letter_code': 'Ser',
        'amino_acid': 'Serine'
    },

    # Threonine
    'ACA': {
        'one_letter_code': 'T',
        'three_letter_code': 'Thr',
        'amino_acid': 'Threonine'
    },
    'ACC': {
        'one_letter_code': 'T',
        'three_letter_code': 'Thr',
        'amino_acid': 'Threonine'
    },
    'ACG': {
        'one_letter_code': 'T',
        'three_letter_code': 'Thr',
        'amino_acid': 'Threonine'
    },
    'ACT': {
        'one_letter_code': 'T',
        'three_letter_code': 'Thr',
        'amino_acid': 'Threonine'
    },

    # Valine
    'GTA': {
        'one_letter_code': 'V',
        'three_letter_code': 'Val',
        'amino_acid': 'Valine'
    },
    'GTC': {
        'one_letter_code': 'V',
        'three_letter_code': 'Val',
        'amino_acid': 'Valine'
    },
    'GTG': {
        'one_letter_code': 'V',
        'three_letter_code': 'Val',
        'amino_acid': 'Valine'
    },
    'GTT': {
        'one_letter_code': 'V',
        'three_letter_code': 'Val',
        'amino_acid': 'Valine'
    },

    # Tryptophan
    'TGG': {
        'one_letter_code': 'W',
        'three_letter_code': 'Trp',
        'amino_acid': 'Tryptophan'
    },

    # Tyrosine
    'TAT': {
        'one_letter_code': 'Y',
        'three_letter_code': 'Tyr',
        'amino_acid': 'Tyrosine'
    },
    'TAC': {
        'one_letter_code': 'Y',
        'three_letter_code': 'Tyr',
        'amino_acid': 'Tyrosine'
    },

    # Stop Codon
    'TAA': {
        'one_letter_code': '*',
        'three_letter_code': 'Ter',
        'amino_acid': 'Stop Codon'
    },
    'TAG': {
        'one_letter_code': '*',
        'three_letter_code': 'Ter',
        'amino_acid': 'Stop Codon'
    },
    'TGA': {
        'one_letter_code': '*',
        'three_letter_code': 'Ter',
        'amino_acid': 'Stop Codon'
    },

    # # Asparagine
    # 'AAC': {
    #     'one_letter_code': 'B',
    #     'three_letter_code': 'Asx',
    #     'amino_acid': ['Asparagine', 'Aspartic Acid']
    # },
    # 'AAT': {
    #     'one_letter_code': 'B',
    #     'three_letter_code': 'Asx',
    #     'amino_acid': ['Asparagine', 'Aspartic Acid']
    # },
    # 'GAC': {
    #     'one_letter_code': 'B',
    #     'three_letter_code': 'Asx',
    #     'amino_acid': ['Asparagine', 'Aspartic Acid']
    # },
    # 'GAT': {
    #     'one_letter_code': 'B',
    #     'three_letter_code': 'Asx',
    #     'amino_acid': ['Asparagine', 'Aspartic Acid']
    # },

    # # Glutamine
    # 'CAA': {
    #     'one_letter_code': 'Z',
    #     'three_letter_code': 'Glx',
    #     'amino_acid': ['Glutamine', 'Glutamic acid']
    # },
    # 'CAG': {
    #     'one_letter_code': 'Z',
    #     'three_letter_code': 'Glx',
    #     'amino_acid': ['Glutamine', 'Glutamic acid']
    # },
    # 'GAA': {
    #     'one_letter_code': 'Z',
    #     'three_letter_code': 'Glx',
    #     'amino_acid': ['Glutamine', 'Glutamic acid']
    # },
    # 'GAG': {
    #     'one_letter_code': 'Z',
    #     'three_letter_code': 'Glx',
    #     'amino_acid': ['Glutamine', 'Glutamic acid']
    # },

    # # X
    # 'NNN': {
    #     'one_letter_code': 'X',
    #     'three_letter_code': 'X',
    #     'amino_acid': ['None']
    # },

}

CODONS = {
    # Alanine
    'GCA': 'Ala',
    'GCC': 'Ala',
    'GCG': 'Ala',
    'GCT': 'Ala',

    # Cysteine
    'TGC': 'Cys',
    'TGT': 'Cys',

    # Aspartic Acid
    'GAC': 'Asp',
    'GAT': 'Asp',

    # Glutamic Acid
    'GAA': 'Glu',
    'GAG': 'Glu',

    # Phenylalanine
    'TTC': 'Phe',
    'TTT': 'Phe',

    # Glycine
    'GGA': 'Gly',
    'GGC': 'Gly',
    'GGG': 'Gly',
    'GGT': 'Gly',

    # Histidine
    'CAC': 'His',
    'CAT': 'His',

    # Isoleucine
    'ATA': 'Ile',
    'ATC': 'Ile',
    'ATT': 'Ile',

    # Lysine
    'AAA': 'Lys',
    'AAG': 'Lys',

    # Leucine
    'CTA': 'Leu',
    'CTC': 'Leu',
    'CTG': 'Leu',
    'CTT': 'Leu',
    'TTA': 'Leu',
    'TTG': 'Leu',

    # Methionine (Start Codon)
    'ATG': 'Met',

    # Asparagine
    'AAC': 'Asn',
    'AAT': 'Asn',

    # Proline
    'CCA': 'Pro',
    'CCC': 'Pro',
    'CCG': 'Pro',
    'CCT': 'Pro',

    # Glutamine
    'CAA': 'Gln',
    'CAG': 'Gln',

    # Arginine
    'AGA': 'Arg',
    'AGG': 'Arg',
    'CGA': 'Arg',
    'CGC': 'Arg',
    'CGG': 'Arg',
    'CGT': 'Arg',

    # Serine
    'AGC': 'Ser',
    'AGT': 'Ser',
    'TCA': 'Ser',
    'TCC': 'Ser',
    'TCG': 'Ser',
    'TCT': 'Ser',

    # Threonine
    'ACA': 'Thr',
    'ACC': 'Thr',
    'ACG': 'Thr',
    'ACT': 'Thr',

    # Valine
    'GTA': 'Val',
    'GTC': 'Val',
    'GTG': 'Val',
    'GTT': 'Val',

    # Tryptophan
    'TGG': 'Trp',

    # Tyrosine
    'TAT': 'Tyr',
    'TAC': 'Tyr',

    # Stop Codon
    'TAA': '*',
    'TAG': '*',
    'TGA': '*'
}


AA_CODES_ONE_TO_THREE = {
    # Alanine
    'A': 'Ala',

    # Cysteine
    'C': 'Cys',

    # Aspartic Acid
    'D': 'Asp',

    # Glutamic Acid
    'E': 'Glu',

    # Phenylalanine
    'F': 'Phe',

    # Glycine
    'G': 'Gly',

    # Histidine
    'H': 'His',

    # Isoleucine
    'I': 'Ile',

    # Lysine
    'K': 'Lys',

    # Leucine
    'L': 'Leu',

    # Methionine (Start Codon)
    'M': 'Met',

    # Asparagine
    'N': 'Asn',

    # Proline
    'P': 'Pro',

    # Glutamine
    'Q': 'Gln',

    # Arginine
    'R': 'Arg',

    # Serine
    'S': 'Ser',

    # Threonine
    'T': 'Thr',

    # Valine
    'V': 'Val',

    # Tryptophan
    'W': 'Trp',

    # Tyrosine
    'Y': 'Tyr',

    # Stop
    'X': 'Ter',
    '*': 'Ter'
}


AMINO_ACIDS = defaultdict(set)
for (codon, aa) in CODONS.items():
    AMINO_ACIDS[aa].add(codon)


def get(codon):
    return CODONS_FULL.get(codon.upper())


def get_codons(three_letter_code):
    return AMINO_ACIDS.get(three_letter_code)


def get_amino_acid(codon):
    return CODONS.get(codon.upper())


# codons that code for the same Amino Acid
def get_isomorphic_codons(codon):
    return get_codons(get_amino_acid(codon)) - set([codon])


def get_stop_codons():
    return set(get_codons('*'))


def codon(_codon):
    return _codon.upper() in CODONS_FULL


def synonymous(codon0, codon1):
    return \
        codon(codon0) and codon(codon1) and \
        get(codon0).get('amino_acid') == get(codon1).get('amino_acid')


def nonsynonymous(codon0, codon1):
    return not synonymous(codon0, codon1)


def stop_codon(_codon):
    return codon(_codon) and _codon.upper() in {'TAA', 'TAG', 'TGA'}


def start_codon(_codon):
    return codon(_codon) and _codon.upper() == 'ATG'


def initiator_codon(_codon):
    return codon(_codon) and _codon.upper() == 'ATG'


def split_into_codons(sequence):
    return [sequence[i:i + 3] for i in range(0, len(sequence), 3)] or ['']
