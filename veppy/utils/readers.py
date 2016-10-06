import csv
import gzip
import logging
from cStringIO import StringIO

logger = logging.getLogger('veppy')


class FileReader(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self._file = None
        self._line_number = -1
        self._logger = logger

    @property
    def file(self):
        if self._file is None:
            is_gzip = self.filepath.endswith('.gz')
            self._file = gzip.open(self.filepath) if is_gzip \
                else open(self.filepath)
            self._logger.debug('reading file \'%s\'.' % self._file.name)
        return self._file

    def __iter__(self):
        return self

    # for use as a context manager
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
        return False

    def __str__(self):
        return "{md}.{cl} (w: {writer})".format(
            md=self.__class__.__module__,
            cl=self.__class__.__name__,
            writer=self._writer
        )

    # noop
    def close(self):
        if self._file:
            self._file.close()

    # hack
    def next(self):
        self._line_number += 1
        return self.file.next().rstrip('\n\r')


class CsvReader(FileReader):
    def __init__(
            self,
            filepath,
            delimiter='\t',
            quotechar='"',
            header=False,
            headercols=None,
            header_processor=lambda x: x.replace(' ', '_').lower(),
            comment='#',
            skiprows=None,
            skipcols=None,
            blacklist=None):
        super(CsvReader, self).__init__(filepath)

        # config
        self.delimiter = delimiter
        self.quotechar = quotechar
        self.header = header
        self.headercols = headercols or []
        self.header_processor = header_processor
        self.comment = comment
        self.skiprows = set(skiprows or [])
        self.skipcols = set(skipcols or [])
        self.blacklist = set(blacklist or [''])
        self._headers_processed = False
        self._headers = None
        self._check_headers()

    def _check_headers(self):
        if not self.header and not self.headercols:
            raise Exception('you must manually specify \'headercols\' \
                if \'header\' is False')

        if self.headercols:
            self._headers_processed = True
            self._headers = list(self.headercols)

    def next(self):
        parsed_row = self._process_next(super(CsvReader, self).next())
        while parsed_row is None:
            parsed_row = self._process_next(super(CsvReader, self).next())
        return parsed_row

    def _process_next(self, line):
        if not line or line.startswith(self.comment):
            return

        if self._line_number in self.skiprows:
            return

        # TODO: make encoding a config option
        line = line.decode('latin1').encode('ascii', 'ignore')
        cols = \
            list(
                csv.reader(
                    StringIO(line),
                    delimiter=self.delimiter,
                    quotechar=self.quotechar
                )
            )[0]
        if not self._headers_processed:
            if self.header:
                self._headers = \
                    self._process_headers(
                        [self.header_processor(c) for c in cols]
                    )
            else:
                self._headers = self.headercols
            self._headers_processed = True

        else:
            return self._process_columns(cols)

    # noop
    def _process_headers(self, headers):
        return headers

    def _process_columns(self, columns):
        # header validation
        if len(columns) != len(self._headers):
            msgs = []
            msg = 'invalid csv reader configuration: ' \
                'data dimensions do not match!'
            msg += ' row has %s columns (skipping %s) but header ' \
                'has %s columns!' % \
                (len(columns), len(self.skipcols), len(self._headers))
            msgs.append(msg)
            msgs.append('row data: %s' % columns)
            msgs.append('headers : %s' % self._headers)
            raise Exception('\n'.join(msgs))

        parsed_row = {}
        for idx, col in enumerate(columns):
            if idx in self.skipcols or \
                    self._headers[idx] in self.skipcols:
                continue

            if col in self.blacklist:
                parsed_row[self._headers[idx]] = None

            else:
                parsed_row[self._headers[idx]] = col

        return parsed_row


class GtfReader(CsvReader):
    def __init__(self, filepath, **kwargs):
        super(GtfReader, self).__init__(
            filepath,
            delimiter='\t',
            headercols=[
                'chromosome',
                'source',
                'feature',
                'start',
                'stop',
                'score',
                'strand',
                'frame',
                'info'
            ],
            blacklist=set(['.', ''])
        )

    # TODO: proper CSV parsing...
    #  'val' typically of the form: 'key "value"''
    #  we should use a proper CSV parsing to handle escaping, but it
    #  appears that all Gtf (k, v) 'attribute' tuples have
    #  singled-valued 'value's
    def _info_to_dict(self, info_str):
        d = {}
        for token in [t.strip() for t in info_str.split('; ')]:
            if token:
                k, v = token.split(' ')
                d[k] = v.strip('"')
        return d

    def _process_columns(self, columns):
        parsed_row = super(GtfReader, self)._process_columns(columns)

        # frame as int
        for k in ['start', 'stop']:
            parsed_row[k] = int(parsed_row[k])

        # frame as int
        parsed_row['frame'] = int(parsed_row['frame'] or 0)

        # convert info to dict...
        parsed_row['info'] = self._info_to_dict(parsed_row['info'])

        return parsed_row


class NcbiMapViewReader(CsvReader):
    def __init__(self, filepath, **kwargs):
        super(NcbiMapViewReader, self).__init__(
            filepath,
            delimiter='\t',
            headercols=[
                'tax_id',
                'chromosome',
                'chr_start',
                'chr_stop',
                'chr_orient',
                'contig',
                'ctg_start',
                'ctg_stop',
                'ctg_orient',
                'feature_name',
                'feature_id',
                'feature_type',
                'group_label',
                'transcript',
                'evidence_code'
            ],
            blacklist=set(['.', ''])
        )

    def _process_columns(self, columns):
        parsed_row = super(NcbiMapViewReader, self)._process_columns(columns)

        # frame as int
        for k in ['chr_start', 'chr_stop', 'ctg_start', 'ctg_stop']:
            parsed_row[k] = int(parsed_row[k])

        return parsed_row


# TODO merge this with other SolveBio VCFReaders
class VcfReader(FileReader):

    def __init__(self, filepath, assembly=None, **kwargs):
        super(VcfReader, self).__init__(filepath)
        self.__reader = None

    # maintain for backwards compatability!
    # TODO: eventually remove and renamed 'self.__reader' to 'self._reader'
    def _reader(self):
        return self.reader

    @property
    def reader(self):
        # TODO this is only used by the test cases
        import vcf
        # add 'strict_whitespace' kwarg to force PyVCF to split
        #  on '\t' only. this has the affect of enabling proper handling
        #  of INFO fields with spaces
        if not self.__reader:
            self.__reader = vcf.Reader(
                filename=self.filepath,
                strict_whitespace=True
            )
        self.build = self.__reader.metadata.get('reference')
        return self.__reader

    @classmethod
    def _get_alternate(cls, obj):
        # if alt is '.' in VCF, PyVCF returns [None]
        #  make that ['.'] instead :-)
        if not obj:
            return '.'

        # if obj has 'sequence' field, it's of class 'vcf.model._Substitution'
        #  which represents a base-case, POJO alternate allele.
        # else, it's a structural variant or something else and
        # we want to skip it
        return getattr(obj, 'sequence', '-')

    def vcf_to_dict(self, vcf_obj):
        return {
            'build': self.build,
            'id': vcf_obj.ID,
            'start': vcf_obj.POS,
            'stop': vcf_obj.POS + len(vcf_obj.REF) - 1,
            'chromosome': vcf_obj.CHROM,
            'reference_allele': vcf_obj.REF,
            'alternate_alleles': map(self._get_alternate, vcf_obj.ALT),
            'info': vcf_obj.INFO
        }

    def __iter__(self):
        return self

    def next(self):
        return self.vcf_to_dict(self.reader.next())
