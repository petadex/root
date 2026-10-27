#!/bin/bash
set -e
source cfg.sh
domain=$1
f=${domain}_pazy_catalytic_domains.fa
f2=${domain}_blastnr_catalytic_domains.fa
f3=${domain}_logan_catalytic_domains.fa
seqtk subseq $data_dir/pazy/$f <(cut -f1 pazy_${domain}.clust | sort -u) > reps/pazy/$f
seqtk subseq $data_dir/nr/$f2 <(cut -f1 nr_${domain}.clust | sort -u) > reps/nr/$f2
seqtk subseq $data_dir/logan/$f3 <(cut -f1 logan_vs_nr_${domain}.unmapped2.clust | sort -u) > reps/logan/$f3
