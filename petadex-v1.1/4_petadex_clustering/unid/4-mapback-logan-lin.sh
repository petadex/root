#!/bin/bash
set -e
source cfg.sh
domain=$1
f=$data_dir/logan/${domain}_logan_catalytic_domains.fa
diamond blastp -d pazy_nr_$domain.reps.faa -q $f --fast -c1 -k1 -o logan_vs_nr_$domain.aln \
	--linsearch --ext full --comp-based-stats 0 --masking 0 --log --algo 0 \
	-f 6 sseqid qseqid --query-cover 90 --id $id_cutoff
cat $f | fasta2tsv.sh | cut -f1 | sort.sh > logan_vs_nr_$domain.ids
comm -2 -3 logan_vs_nr_$domain.ids <(cut -f2 logan_vs_nr_${domain}.aln | sort.sh) > logan_vs_nr_$domain.unmapped
seqtk subseq $f logan_vs_nr_$domain.unmapped > logan_vs_nr_$domain.unmapped.faa
rm logan_vs_nr_$domain.ids
rm logan_vs_nr_$domain.unmapped
