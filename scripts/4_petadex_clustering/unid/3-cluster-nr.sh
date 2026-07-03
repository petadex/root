#!/bin/bash
source cfg.sh
f="nr_vs_pazy_$1.unmapped.faa"
/root/tmp/clustering/unid/cluster-full.sh $1 nr $f
cat pazy_$1.reps.faa nr_$1.reps.faa > pazy_nr_$1.reps.faa
