import os
import logging

from collections import defaultdict
from operator import attrgetter

from veppy import SUPPORTED_BUILDS
from veppy import DATA_DIR
from veppy import SOURCE_DATA

from .data import CGD_GENE_SYMBOLS

logger = logging.getLogger('veppy')


# retun genes and transcripts in "natural" sort order
def sort_genes(genes):
    return sorted(genes, key=lambda g: g.length, reverse=True)


# NOTE: assumes all transcripts are within the same gene
def sort_transcripts(transcripts):
    # yields  transcript dict sorted by length and then by accession number
    # this is logic is very important to maintain!
    return sorted(
        sorted(
            transcripts,
            key=attrgetter('accession_version'), reverse=True
        ),
        key=attrgetter('accession_int')
    )


def get_default_gene(genes):
    genes_in_cgd = [g for g in genes if g.symbol in CGD_GENE_SYMBOLS]

    if len(genes_in_cgd) == 1:
        default_gene = genes_in_cgd[0]

    elif len(genes_in_cgd) > 1:
        default_gene = sort_genes(genes_in_cgd)[0]
        logger.info(
            'Found multiple genes %s in CGD. '
            'Fallback to using the longest gene %s.' %
            ([g.symbol for g in genes_in_cgd], default_gene.symbol)
        )

    else:
        # len(genes_in_cgd) can be 0 because either no CGD genes were laoded
        #  OR because there simply was no overlap between
        #  'genes' and 'CGD_GENE_SYMBOLS'
        default_gene = sort_genes(genes)[0]
        logger.info(
            'Could not find genes %s in CGD. '
            'Fallback to using the longest gene %s. ' %
            ([g.symbol for g in genes], default_gene.symbol)
        )

    return default_gene


# assumes all transcripts are from the SAME gene
def get_default_transcript(transcripts):
    return sort_transcripts(transcripts)[0]


def set_defaults(transcripts):
    # group transcripts by Gene and find defaults
    _groups = defaultdict(list)
    for transcript in transcripts:
        _groups[transcript.gene].append(transcript)

    default_gene = get_default_gene(_groups.keys())
    if default_gene:
        default_gene.default = True

    for (gene, transcripts) in _groups.items():
        default_transcript = get_default_transcript(transcripts)
        if default_transcript:
            default_transcript.default = True


def build_data_directory(build):
    release, version = parse_build(build)
    return os.path.join(DATA_DIR, release)


def fasta_filepath(build, source):
    return os.path.join(
        build_data_directory(build),
        SOURCE_DATA[build]['sequence'][source]
    )


def feature_filepath(build, gene_model):
    return os.path.join(
        build_data_directory(build),
        SOURCE_DATA[build]['feature'][gene_model]
    )


def parse_build(build):
    if '.' in build:
        return build.split('.', 1)
    else:
        return (build, '')


def supported_builds():
    return list(SUPPORTED_BUILDS)


def supported_build(build):
    release, version = parse_build(build)
    for b in SUPPORTED_BUILDS:
        if b.upper() == release.upper():
            return b

    logger.error('Unsupported genome build %s.' % build)
    return None
