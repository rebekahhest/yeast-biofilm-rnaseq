#!/bin/bash

# Download FASTQ files
mkdir -p data/raw
cd data/raw
for i in {57..65}; 
do
  echo "Downloading SRR105516${i}"
  prefetch SRR105516${i}
  fasterq-dump SRR105516${i}
  rm -r SRR105516${i}
  
done
cd ../../