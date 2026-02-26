#!/bin/bash

# Download FASTQ files
mkdir -p data/raw
cd data/raw
for i in {57..65}; 
do
  echo "Downloading SRR105516${i}"
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/0${i}/SRR105516${i}/SRR105516${i}.fastq.gz  
done
cd ../../