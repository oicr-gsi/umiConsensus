#!/bin/bash
set -euo pipefail

cd "$1"
module load samtools/1.9 vcftools/0.1.16 1>/dev/null
ls | sed 's/.*\.//' | sort | uniq -c

find . -name "*.stats.txt" -xtype f -exec sh -c '
    f="{}"
    first_line=$(head -n1 "$f")

    if echo "$first_line" | grep -qE "^[A-Za-z]"; then
        # First line looks like a header (starts with a letter) —
        # keep it pinned, sort the rest
        (echo "$first_line"; tail -n +2 "$f" | LC_ALL=C sort)
    else
        # No header detected — sort the whole file
        LC_ALL=C sort "$f"
    fi | md5sum
' \;

find . -name "*.bam" -xtype f -exec sh -c "samtools flagstat {} | grep 'in total' | md5sum | sort" \;