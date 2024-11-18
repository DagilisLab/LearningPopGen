#!/bin/bash

gene=$1
taxon=$2

Rscript download_data.R -gene=$gene -taxon=$taxon

mafft --auto data/${taxon}_${gene}.fasta > data/aligned_${taxon}_${gene}.fasta

iqtree2 -s data/aligned_${taxon}_${gene}.fasta 


