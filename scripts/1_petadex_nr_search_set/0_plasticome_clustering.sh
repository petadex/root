usearch --cluster_fast plasticome.preclust.fa \
  -id 0.30 \
  -centroids plasticome.30.fa \
  -uc plasticome.30.uc

# usearch v11.0.667_i86linux64, 65.5Gb RAM, 20 cores
# (C) Copyright 2013-18 Robert C. Edgar, all rights reserved.
# https://drive5.com/usearch

# License: a.babaian@utoronto.ca, non-profit use, max 1 process(es)

# 00:00 43Mb    100.0% Reading plasticome.preclust.fa
# 00:00 9.9Mb  CPU has 20 cores, defaulting to 10 threads
# 00:00 547Mb   100.0% DF
# 00:00 689Mb  213 seqs, 205 uniques, 198 singletons (96.6%)
# 00:00 689Mb  Min size 1, median 1, max 3, avg 1.04
# 00:00 695Mb   100.0% DB
# 00:00 762Mb   100.0% 89 clusters, max size 24, avg 2.4
# 00:00 762Mb   100.0% Writing centroids to plasticome.30.fa
                                                          
#       Seqs  205
#   Clusters  89
#   Max size  24
#   Avg size  2.4
#   Min size  1
# Singletons  61, 29.8% of seqs, 68.5% of clusters
#    Max mem  762Mb
#       Time  1.00s
# Throughput  205.0 seqs/sec.

# All vs. All Pairwise Alignment
#query target id alnlen mism opens qlo qhi tlo thi evalue bits
usearch -allpairs_local plasticome.preclust.fa \
        -id 0.25 -evalue 0.1 \
        -blast6out plasticome.aln


## Install BLAST
# wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz ./
# tar -xvf ncbi-blast-2.16.0+-x64-linux.tar.gz
#
# cp ncbi-blast-2.16.0+/bin/* ~/bio/bin/
## Download BLAST nr
# mkdir -p nr; cd nr
# update_blastdb.pl --source aws --num_threads 4 nr
# diamond prepdb -d nr

# Run diamond blastp
~/bio/bin/diamond blastp \
  -q plasticome.v1.fa \
  -d /data/nr/nr \
  --very-sensitive \
  -p16 -k0 \
  -f 6 qseqid  qstart qend qlen \
       sseqid  sstart send slen \
       pident evalue bitscore \
  > plasticome.blast.fa

# Retrieve sequence hits