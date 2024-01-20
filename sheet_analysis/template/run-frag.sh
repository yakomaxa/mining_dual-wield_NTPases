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


# link the original foldcomp database to the splitdb. This contains full data (for P-loop containing subset).
ln /lab/work/ksakuma/myprojects/FOLDCOMP/parallel_UniProt/gen02/UniProtV4WalkerLike/uniprotv4walkerlike ./splitdb

# link the chuncked index file given as ${1} to splitdb. The name has to be uniprotv4walkerlike.index to deceive foldcomp.
# This let foldcomp to only read the part of input db acoording to the chunked index.
ln ${1} ./splitdb/uniprotv4walkerlike.index

# Run my script
/lab/work/ksakuma/myprojects/FOLDCOMP/venv/bin/python ./afdb_nishi.py
