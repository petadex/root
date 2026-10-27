#!/bin/bash
set -e
domain=$1
merged=logan_vs_nr_$domain.unmapped2.clust_merged
merge.sh logan_vs_nr_$domain.unmapped2.clust logan_vs_nr_$domain.unmapped2.clust2 > $merged
