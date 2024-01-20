#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -pe OpenMP 4
export OMP_NUM_THREADS=1

mkdir ./structures
mkdir ./splitdb
mkdir ./fastas
mkdir ./fragments
mkdir ./fragment_fastas

# link the original foldcomp database to the splitdb. This contains full data.
ln /lab/work/ksakuma/myprojects/FOLDCOMP/UniProt/afdb_uniprot_v4 ./splitdb

# link the chuncked index file given as ${1} to splitdb. The name has to be afdb_uniprot_v4.index to deceive foldcomp.
# This let foldcomp to only read the part of db in the chunked index.
ln ${1} ./splitdb/afdb_uniprot_v4.index

# Run my script
/lab/work/ksakuma/myprojects/FOLDCOMP/venv/bin/python ./afdb_all_4.py
