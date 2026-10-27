#!/bin/bash
set -e
domain=$1
cat split-$domain/out_* > logan_vs_nr_$domain.unmapped2.self_aln
