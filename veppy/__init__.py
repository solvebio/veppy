import os
import sys
import logging

log_level = logging.INFO
logger = logging.getLogger('veppy')

if not logger.handlers:
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    formatter = logging.Formatter(
        '%(asctime)s - [%(name)s (%(process)d)] - %(levelname)s: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

# Veppy Settings
DATA_DIR = os.environ.get(
    'VEPPY_DATA_DIR',
    './data'
    # os.path.abspath(
    #     os.path.join(os.path.dirname(__file__), '..', 'data')
    # )
)

SOURCE_DATA = {
    'GRCh37': {
        'sequence': {
            'genbank': 'Homo_sapiens.GRCh37.75.genbank.fa'
        },
        'feature': {
            # refseq == NCBI alias for now until we transition to using refseq.
            'refseq': 'seq_gene.104.105.combined.sorted.md',
            'ncbi': 'seq_gene.104.105.combined.sorted.md',
            'gencode': 'gencode.v19.annotation.gtf'
        }
    }
}

SOURCE_DATA_MD5 = {
    'Homo_sapiens.GRCh37.75.genbank.fa':
    '1a4c08a7a4be4cce8bdfe06f71802319',
    'Homo_sapiens.GRCh37.75.genbank.fa.flat':
    '3821935a6f66fd7fbe5a150701762ed0',
    'Homo_sapiens.GRCh37.75.genbank.fa.gdx':
    '824390ca5a69f788208723008f7dbfca'
}

# fail if md5sum fails
MD5_CHECK_FAIL = os.environ.get('MD5_CHECK_FAIL', False)

SUPPORTED_BUILDS = os.environ.get(
    'SUPPORTED_BUILDS', 'GRCh37').split(',')

SUPPORTED_GENE_MODELS = SOURCE_DATA['GRCh37']['feature'].keys()

FEATURE_LIMIT = int(os.getenv('FEATURE_LIMIT', sys.maxint))


# Testing

# use for skipping tests that require full data
RUN_FULL_DATA_TESTS = \
    os.environ.get('RUN_FULL_DATA_TESTS', True) is True
