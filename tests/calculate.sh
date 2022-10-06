#!/bin/bash
set -euo pipefail

cd $1
module load samtools/1.9 vcftools/0.1.16 1>/dev/null
ls | sed 's/.*\.//' | sort | uniq -c
find . -name "*.stats.txt" -xtype f -exec sh -c "cat {} | md5sum" \;
find . -name "*.bam" -xtype f -exec sh -c "samtools flagstat {} | grep 'in total' | md5sum | sort" \;
