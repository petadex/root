#!/bin/bash
set -e
source cfg.sh
domain=$1
f=logan_vs_nr_$domain.unmapped2.faa
clust=logan_vs_nr_$domain.unmapped2.clust
diamond deepclust -d $f -o $clust --cluster-steps $steps_logan_linclust -M 2000G \
        --ext full --comp-based-stats 0 --masking 0 --log --id $id_cutoff --member-cover 90

#seqtk subseq $f <(cut -f1 $clust | uniq) > logan_vs_nr_$domain.unmapped2.clust.reps.faa
