#!/bin/bash
source cfg.sh
domain=$1
f="$data_dir/pazy/${domain}_pazy_catalytic_domains.fa"
echo $f
echo $id_cutoff
/root/tmp/clustering/unid/cluster-all-vs-all.sh $domain pazy $f
