#!/bin/bash
DATA_DIR='data/GRCh37'
mkdir -p $DATA_DIR

echo "Downloading reference data to ${DATA_DIR}"

GZFILE="gencode.v19.annotation.gtf.gz"
OUTFILE="gencode.v19.annotation.gtf"
if [ -e "$DATA_DIR/$OUTFILE" ]; then
    echo "Skipping ${GZFILE} since $DATA_DIR/$OUTFILE exists..."
else
    echo "Downloading ${GZFILE} ..."
    curl -L "http://veppy.s3-website-us-east-1.amazonaws.com/${GZFILE}" > ${DATA_DIR}/${GZFILE}
    gunzip ${DATA_DIR}/${GZFILE}
fi

GZFILE="seq_gene.104.105.combined.sorted.md.gz"
OUTFILE="seq_gene.104.105.combined.sorted.md"
if [ -e "$DATA_DIR/$OUTFILE" ]; then
    echo "Skipping ${GZFILE} since $DATA_DIR/$OUTFILE exists..."
else
    echo "Downloading ${GZFILE} ..."
    curl -L "http://veppy.s3-website-us-east-1.amazonaws.com/${GZFILE}" > ${DATA_DIR}/${GZFILE}
    gunzip ${DATA_DIR}/${GZFILE}
fi

OUTFILE='Homo_sapiens.GRCh37.75.genbank.fa'
if [ -e "$DATA_DIR/$OUTFILE" ]; then
    echo "Skipping ${GZFILE} since $DATA_DIR/$OUTFILE exists..."
else
    CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

    for i in $CHROMS
    do
        echo "Downloading chromosome ${i}";
        curl http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/chr${i}.fa.gz | gzip -d >> ${DATA_DIR}/${OUTFILE}
    done
fi
