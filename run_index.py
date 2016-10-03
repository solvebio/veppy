from veppy.utils.fasta import GenbankFastaFile
from veppy.feature_file import GencodeGtfFile
from veppy.feature_file import NcbiMapViewFile

print "Indexing FASTA and feature files, this may take 10 minutes or more..."

fasta_file = GenbankFastaFile.for_build("GRCh37")
feature_file_gencode = GencodeGtfFile.for_build('GRCh37')
feature_file_NCBI = NcbiMapViewFile.for_build('GRCh37')
