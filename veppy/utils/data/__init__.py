import os
import gzip
import pickle

data_dir = os.path.dirname(__file__)

CGD_GENE_SYMBOLS_VERSION = "CGD/1.1.2-2015-09-24/CGD"
with gzip.open(os.path.join(
        data_dir, 'cgd-1.1.2-2015-09-24-cgd.pickle.gz'), 'rb') as f:
    CGD_GENE_SYMBOLS = pickle.load(f)
