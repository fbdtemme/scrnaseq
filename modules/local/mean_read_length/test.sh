#!/usr/bin/env bash
reads=()
reads+=( S10_L001_R2_001.fastq.gz SRR8599150_S1_L001_R2_001.sub5000.fastq.gz )

count=0;
sum=0;

for r in $reads; do
    echo $r
    r=\$( seqkit stats -T "\$r" |  sed -n 2p | cut -d\$'\t' -f7 | xargs echo -n );
    (( sum+=r ))
    (( ++count ));
done;

echo \$(( sum / count ))