#!/bin/bash
source cfg.sh
f="$data_dir/nr/$1_blastnr_catalytic_domains.fa"
diamond blastp -q $f -d pazy_$1.reps.faa -c1 -k1 -o nr_vs_pazy_$1.aln \
	--ext full --comp-based-stats 0 --masking 0 --log --algo 0 \
	-f 6 sseqid qseqid --id $id_cutoff --query-cover 90 --ultra-sensitive
cat $f | fasta2tsv.sh | cut -f1 | sort.sh > nr_vs_pazy_$1.ids
comm -2 -3 nr_vs_pazy_$1.ids <(cut -f2 nr_vs_pazy_$1.aln | sort.sh) > nr_vs_pazy_$1.unmapped
seqtk subseq $f nr_vs_pazy_$1.unmapped > nr_vs_pazy_$1.unmapped.faa
rm nr_vs_pazy_$1.ids
rm nr_vs_pazy_$1.unmapped
