#! /usr/bin/env bash

set -e

# hard-parallelize ukbb ldsc jobs
# - loops specified number of batches
# - spins up cluster
# - submits ldsc batch in background

maxi=$((3))

for i in `seq 1 $maxi`; do
# for i in `seq 1 1`; do
	
	cluster start ukbb-rkw${i} -m n1-standard-16 --num-workers 0 --num-preemptible-workers 0
	cluster submit ukbb-rkw${i} ldsc_h2part_parallel_batch.py --args "--parsplit ${maxi} --paridx ${i}" &

done


# eof
