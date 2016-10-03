import os
import logging

from abc import ABCMeta, abstractmethod

from pyfasta import Fasta, FastaRecord
from pyfasta.fasta import FastaNotFound

from .helpers import fasta_filepath
from .helpers import supported_build
from .helpers import parse_build
from .md5 import calculate_md5

from .. import SUPPORTED_BUILDS
from .. import SOURCE_DATA_MD5
from .. import MD5_CHECK_FAIL
from .. import errors

logger = logging.getLogger('veppy')


class FastaFile(object):
    __metaclass__ = ABCMeta

    # -- Util methods
    @classmethod
    def load_build(klass, build):
        try:
            translated_build = supported_build(build)
            if not klass.BUILD_CACHE.get(translated_build):
                klass.BUILD_CACHE[translated_build] = \
                    klass(build, klass.filepath_for_build(translated_build))
        except (FastaNotFound, errors.FastaFileException), e:
            logger.fatal(
                'Error loading FastA file for build %s: %s' % (build, e)
            )
            raise

    @classmethod
    def for_build(klass, build):
        _build = supported_build(build)
        if not _build:
            raise errors.VeppyFileException(
                'Unsupported build \'%s\'.' % build)

        klass.load_build(_build)
        return klass.BUILD_CACHE[_build]

    @classmethod
    @abstractmethod
    def filepath_for_build(klass, build):
        pass

    def index(self, reindex=False):
        return self.get_sequence('1', 1, 1)
    # ---

    def __init__(self, build, filepath, index=False):
        self.build = build
        self.filepath = filepath

        # check for FastA file...
        if not os.path.exists(self.filepath):
            raise errors.FastaFileException(
                'FastA file \'%s\' not found!' % self.filepath
            )

        # check for FastA file... index files FIRST (!!) because Fasta(...)
        #  __init__ checks for indexes and greedily creates
        #  them if they don't exist...
        if not os.path.exists(self.filepath + '.flat') or \
                not os.path.exists(self.filepath + '.gdx'):
            logger.warn(
                'FASTA indexes not found. Indexing...\'%s\'. '
                'This may take 10 minutes or more.' % self.filepath
            )

        self._fasta = Fasta(self.filepath, record_class=FastaRecord,
                            flatten_inplace=True)

        # Ensure the MD5s are correct
        # TODO: This is slow (3s per file). Move somewhere else.
        files = [self.filepath,
                 self.filepath + '.flat',
                 self.filepath + '.gdx']
        for f in files:
            calculated = calculate_md5(f)
            expected = SOURCE_DATA_MD5.get(os.path.basename(f))
            if calculated != expected:
                message = \
                    'FASTA file {} did not pass MD5 check. ' \
                    'Expected: {} Found: {}' \
                    .format(f, expected, calculated)
                if MD5_CHECK_FAIL:
                    raise errors.FastaFileException(message)
                else:
                    logger.warn(message)
            else:
                logger.info('Fasta file {} passed MD5 check.'
                            .format(f))

        # TODO: eventually, this should be dynamic and file-specific
        self.chromosomes = \
            map(lambda x: str(x), range(1, 23)) + ['MT', 'X', 'Y']

    def get(self, chromosome, start, stop):
        fasta_chromosome = self.get_chromosome(chromosome)
        if not fasta_chromosome:
            logger.warning('Invalid chromosome requested in FASTA reader: {}'
                           .format(chromosome))
            return None

        # protect pyfasta from itself...if start > stop, PyFasta will
        #  happily fetch the wrapping range of [stop:] + [:start]
        #  and more than likely OOM you into oblivion in the process.
        if start > stop:
            logger.debug(
                'Invalid range found. chr: %s, range: [%s, %s].' %
                (chromosome, start, stop)
            )
            return ''

        return fasta_chromosome[start - 1:stop - 1 + 1]

    # backwards compat w/ Sequence API
    def get_sequence(self, chromosome, start, stop):
        return self.get(chromosome, start, stop)

    def get_slice(self, chromosome, _slice):
        return self.get(chromosome, _slice.start, _slice.stop - 1)

    @abstractmethod
    def get_chromosome(self, chromosome):
        pass

    @abstractmethod
    def get_chromosome_accession(self, chromosome):
        pass

    @staticmethod
    def supported_chromosome_accession(chromosome_accession):
        return chromosome_accession.upper() in \
            REFSEQ_CHROMOSOME_ACCESSIONS

    @staticmethod
    def get_build_from_chromosome_accession(chromosome_accession):
        if FastaFile.supported_chromosome_accession(chromosome_accession):
            return REFSEQ_CHROMOSOME_ACCESSIONS[
                chromosome_accession.upper()]['build']
        return None

    @staticmethod
    def get_chromosome_from_chromosome_accession(chromosome_accession):
        if FastaFile.supported_chromosome_accession(chromosome_accession):
            return REFSEQ_CHROMOSOME_ACCESSIONS[
                chromosome_accession.upper()]['chromosome']
        return None


# Genbank sequence Fasta file
class GenbankFastaFile(FastaFile):
    BUILD_CACHE = {k: None for k in SUPPORTED_BUILDS}

    FASTA_CHROMOSOME_ACCESSIONS = {
        'GRCH37': {
            'CM000663.1': 'NC_000001.10',
            'CM000664.1': 'NC_000002.11',
            'CM000665.1': 'NC_000003.11',
            'CM000666.1': 'NC_000004.11',
            'CM000667.1': 'NC_000005.9',
            'CM000668.1': 'NC_000006.11',
            'CM000669.1': 'NC_000007.13',
            'CM000670.1': 'NC_000008.10',
            'CM000671.1': 'NC_000009.11',
            'CM000672.1': 'NC_000010.10',
            'CM000673.1': 'NC_000011.9',
            'CM000674.1': 'NC_000012.11',
            'CM000675.1': 'NC_000013.10',
            'CM000676.1': 'NC_000014.8',
            'CM000677.1': 'NC_000015.9',
            'CM000678.1': 'NC_000016.9',
            'CM000679.1': 'NC_000017.10',
            'CM000680.1': 'NC_000018.9',
            'CM000681.1': 'NC_000019.9',
            'CM000682.1': 'NC_000020.10',
            'CM000683.1': 'NC_000021.8',
            'CM000684.1': 'NC_000022.10',
            'CM000685.1': 'NC_000023.10',
            'CM000686.1': 'NC_000024.9'
        }
    }

    @classmethod
    def filepath_for_build(klass, build):
        (release, version) = parse_build(build)
        logger.info(
            'Loading build %s (release: %s, version: %s).' %
            (build, release, version or 'None')
        )

        return fasta_filepath(build, 'genbank')

    # chromosomes entries in *.fa files look like:
    #    'gi|224384768|gb|CM000663.1| Homo sapiens chromosome 1, GRC primary reference assembly'            # noqa
    # given an input chromosome 'chrN', strip leading 'ch(r)' and lookup based
    #  on keys that match 'N '
    def _get_fasta_key_for_chromosome(self, chromosome):
        chr_clean = chromosome.upper().lstrip('CHR')
        for k in self._fasta.keys():
            if k.split(', ')[0].upper().endswith(' %s' % chr_clean):
                return k
        return None

    def get_chromosome(self, chromosome):
        k = self._get_fasta_key_for_chromosome(chromosome)
        return self._fasta.get(k) if k else None

    def get_chromosome_accessions(self):
        return self.REFSEQ_CHROMOSOME_ACCESSIONS[self.build.upper()]

    def get_chromosome_accession(self, chromosome):
        fasta_key = self._get_fasta_key_for_chromosome(chromosome)
        if fasta_key:
            chr_accession = fasta_key.split('|')[3]
            return self.FASTA_CHROMOSOME_ACCESSIONS[self.build.upper()] \
                .get(chr_accession)
        return None


REFSEQ_CHROMOSOME_ACCESSIONS = {
    'NC_000001.10': {
        'chromosome': '1',
        'build': 'GRCh37'
    },
    'NC_000002.11': {
        'chromosome': '2',
        'build': 'GRCh37'
    },
    'NC_000003.11': {
        'chromosome': '3',
        'build': 'GRCh37'
    },
    'NC_000004.11': {
        'chromosome': '4',
        'build': 'GRCh37'
    },
    'NC_000005.9': {
        'chromosome': '5',
        'build': 'GRCh37'
    },
    'NC_000006.11': {
        'chromosome': '6',
        'build': 'GRCh37'
    },
    'NC_000007.13': {
        'chromosome': '7',
        'build': 'GRCh37'
    },
    'NC_000008.10': {
        'chromosome': '8',
        'build': 'GRCh37'
    },
    'NC_000009.11': {
        'chromosome': '9',
        'build': 'GRCh37'
    },
    'NC_000010.10': {
        'chromosome': '10',
        'build': 'GRCh37'
    },
    'NC_000011.9': {
        'chromosome': '11',
        'build': 'GRCh37'
    },
    'NC_000012.11': {
        'chromosome': '12',
        'build': 'GRCh37'
    },
    'NC_000013.10': {
        'chromosome': '13',
        'build': 'GRCh37'
    },
    'NC_000014.8': {
        'chromosome': '14',
        'build': 'GRCh37'
    },
    'NC_000015.9': {
        'chromosome': '15',
        'build': 'GRCh37'
    },
    'NC_000016.9': {
        'chromosome': '16',
        'build': 'GRCh37'
    },
    'NC_000017.10': {
        'chromosome': '17',
        'build': 'GRCh37'
    },
    'NC_000018.9': {
        'chromosome': '18',
        'build': 'GRCh37'
    },
    'NC_000019.9': {
        'chromosome': '19',
        'build': 'GRCh37'
    },
    'NC_000020.10': {
        'chromosome': '20',
        'build': 'GRCh37'
    },
    'NC_000021.8': {
        'chromosome': '21',
        'build': 'GRCh37'
    },
    'NC_000022.10': {
        'chromosome': '22',
        'build': 'GRCh37'
    },
    'NC_000023.10': {
        'chromosome': 'X',
        'build': 'GRCh37'
    },
    'NC_000024.9': {
        'chromosome': 'Y',
        'build': 'GRCh37'
    }
}
