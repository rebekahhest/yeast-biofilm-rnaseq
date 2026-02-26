#!/bin/bash

# Assess Quality using FASTQC
Mkdir -p data/fastqc

for i in data/raw/SRR105516{57..65}.fastq.gz
do
  samp=$(basename ${i} .fastq.gz)
  echo "Running FASTQC on ${samp}" 
  fastqc ${i} -o data/fastqc
done