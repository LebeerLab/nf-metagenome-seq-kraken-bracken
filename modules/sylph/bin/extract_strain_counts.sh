#!/bin/bash

folder=$1

for f in $folder/*.sylphmpa; do
  sample=${f%%.sylphmpa}
  sample=${sample##./}
  echo $sample
  cat $f | 
  grep t__ | 
  cut -f 1,2 |
  sed -s "s/\s/_/" | 
  sed -e "s/$/\t$sample/g"  >> counts
done

datamash --sort --whitespace --filler=0 crosstab 1,3 sum 2 < counts > count_matrix.tsv
