#!/bin/bash
set -e
source cfg.sh
domain=$1
cat split-$domain/out_* > logan_vs_nr_$domain.aln2
f="logan_vs_nr_$domain.unmapped.faa"
cat $f | fasta2tsv.sh | cut -f1 | sort.sh > logan_vs_nr_$domain.ids2
comm -2 -3 logan_vs_nr_$domain.ids2 <(cut -f2 logan_vs_nr_${domain}.aln2 | sort.sh) > logan_vs_nr_$domain.unmapped2
seqtk subseq $f logan_vs_nr_$domain.unmapped2 > logan_vs_nr_$domain.unmapped2.faa
rm logan_vs_nr_$domain.ids2
rm logan_vs_nr_$domain.unmapped2
