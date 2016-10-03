# veppy (alpha)

**WARNING: This code is an alpha release and not production-ready. APIs may change at any time.**

A genetic variant effect predictor for Python. Inspired by SnpEff and VEP.


## Installation

```
python setup.py install
```

## Setup

### Prepare a directory for veppy data

The default data path is: `./data`

(OPTIONAL) You can override this by setting `$VEPPY_DATA_DIR`.

```
export VEPPY_DATA_DIR=/opt/veppy
```

### Download source data and build indexes

**NOTE: This step downloads about 1gb of data. After indexing, the data directory will consume about 8gb of disk space.**

```
./scripts/download_data_GRCh37.sh
```

### Build indexes (this only needs to be done once)

```
python ./run_index.py
```

## Example Usage

```
from veppy.veppy import calculate_consequences
variant = ('1', 8025384, 'A', 'T')
result = calculate_consequences('GRCh37', *variant)
print result.results
```


## Testing

Tests are based on chr1 versions of input data (TODO upload this somewhere)

```
$ nosetests
```

Coverage:
```
$ nosetests --with-coverage --cover-package=veppy
```


## About SolveBio

SolveBio is a genomics company based in New York City.

[![SolveBio](https://s3.amazonaws.com/veppy/solvebio_logo.png)](https://www.solvebio.com/)
