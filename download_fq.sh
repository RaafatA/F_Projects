#!/bin/bash
# Download paired-end FASTQ files from SRA using prefetch and fasterq-dump
# Each accession file gets its own folder
# Usage: ./download_fastq.sh accession_file1.txt accession_file2.txt ...

if [ $# -lt 1 ]; then
    echo "Usage: $0 accession_file1.txt accession_file2.txt ..."
    exit 1
fi

for ACCESSION_LIST in "$@"; do
    BASENAME=$(basename "$ACCESSION_LIST" .txt)

    echo ">>> Processing file: $ACCESSION_LIST"

    # Make a folder for this accession list
    mkdir -p "$BASENAME/sra_files" "$BASENAME/fastq_files"

    while read -r ACC; do
        # Skip empty lines
        if [[ -z "$ACC" ]]; then
            continue
        fi

        echo "=== Downloading $ACC from $ACCESSION_LIST ==="

        # Step 1: Download SRA file with prefetch
        prefetch -O "$BASENAME/sra_files" "$ACC"

        # Step 2: Convert to FASTQ (paired-end, split into _1 and _2)
        fasterq-dump "$BASENAME/sra_files/$ACC" \
            -O "$BASENAME/fastq_files" \
            --split-files \
            --progress

        echo "Finished $ACC"
    done < "$ACCESSION_LIST"

    echo ">>> Completed file: $ACCESSION_LIST (results in $BASENAME/)"
done

echo "All downloads completed."
