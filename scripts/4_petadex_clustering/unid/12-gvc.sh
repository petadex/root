#!/bin/bash
set -e
domain=$1
input=logan_vs_nr_$domain.unmapped2.clust.reps.faa
aln=logan_vs_nr_$domain.unmapped2.self_aln
clust=logan_vs_nr_$domain.unmapped2.clust2
cat $input | fasta2tsv.sh | cut -f1 > $domain.ids
diamond greedy-vertex-cover --edges $aln -d $domain.ids -o $clust --member-cover 90
