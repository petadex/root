# self alignment unmapped logan reps after linclust

set -e
i=$2
domain=$(sed "$1!d" ../domains)
f=split-$domain/seqs_$i.fasta
echo "Domain=$domain"

diamond blastp -d logan_vs_nr_$domain.unmapped2.clust.reps.faa -q $f --fast -c1 -k0 -o split-$domain/out_$i -b0.2 \
        --ext full --comp-based-stats 0 --masking 0 --algo 0 -t /dev/shm \
        -f 6 qseqid sseqid qcovhsp scovhsp corrected_bitscore --log --query-or-subject-cover 90 --id 90 --id2 17
echo $i >> fin
