#!/bin/bash

# Perform quantification using Salmon
mkdir -p data/yeast_index
mkdir -p data/quant 

# Build index
salmon index \
  -t data/reference/GCF_000146045.2_R64_rna.fna \
  -i data/yeast_index

# Estimate transcript abundances
for i in data/raw/SRR105516{57..65}.fastq.gz;
do
  samp=$(basename ${i} .fastq.gz)
  echo "Quantifying ${samp}"
salmon quant \
  -i data/yeast_index \
  -l A \
  -r ${i} \
  -p 4 \
  --validateMappings \
  -o data/quant/${samp}
done 



