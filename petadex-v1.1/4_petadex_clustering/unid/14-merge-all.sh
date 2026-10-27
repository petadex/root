#!/bin/bash
set -e
source cfg.sh
domain=$1
merged=logan_vs_nr_$domain.unmapped2.clust
cat pazy_$domain.clust nr_vs_pazy_$domain.aln nr_$domain.clust logan_vs_nr_$domain.aln logan_vs_nr_$domain.aln2 $merged | sort.sh -k1,1 > $domain.all.tsv
