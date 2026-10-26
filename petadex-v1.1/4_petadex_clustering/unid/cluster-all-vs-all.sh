#!/bin/bash
set -e
domain=$1
code=$2
file=$3
id=${code}_${domain}
source cfg.sh
echo id_cutoff=$id_cutoff
diamond blastp -q $file -d $file -o $id.aln -f 6 qseqid sseqid qcovhsp scovhsp corrected_bitscore \
	--ext full --comp-based-stats 0 --masking 0 --log --id $id_cutoff --query-or-subject-cover 90 -c1 -k0 --ultra-sensitive
cat $file | fasta2tsv.sh | cut -f1 > $id.ids
diamond greedy-vertex-cover --edges $id.aln -d $id.ids -o $id.clust.tmp --member-cover 90
sort.sh -k1,1 $id.clust.tmp > $id.clust
rm $id.clust.tmp
seqtk subseq $file <(cut -f1 $id.clust | uniq) > $id.reps.faa
rm $id.ids
rm $id.aln
