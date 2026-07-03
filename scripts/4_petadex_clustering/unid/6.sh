#mapback-logan

set -e
source cfg.sh
i=$2
domain=$(sed "$1!d" ../domains)
f=split-$domain/seqs_$i.fasta
echo "Domain=$domain"

diamond blastp -d pazy_nr_$domain.reps.faa -q $f $sens_backmap_logan -c1 -k1 -o split-$domain/out_$i -b1 \
        --ext full --comp-based-stats 0 --masking 0 --algo 0 -t /dev/shm \
        -f 6 sseqid qseqid --log --query-cover 90 --id $id_cutoff
echo $i >> fin