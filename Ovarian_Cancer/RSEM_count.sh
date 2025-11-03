#!/bin/bash

REF="/home/2375894/radwa-scratch/GEO_GTEx/refrence"
OUTDIR="/home/2375894/radwa-scratch/GEO_GTEx/rsem_results"
THREADS=8

# Enable recursive globbing
shopt -s nullglob

for STUDY in /home/2375894/radwa-scratch/GEO_GTEx/mapped_GSE*/; do
  STUDYNAME=$(basename "$STUDY")
  mkdir -p "$OUTDIR/$STUDYNAME"

  echo "========== Processing Study: $STUDYNAME =========="

  # Match STAR coordinate-sorted BAMs
  for BAM in "$STUDY"/*_Aligned.sortedByCoord.out.bam; do
    [ -e "$BAM" ] || { echo "No BAM files in $STUDYNAME"; continue; }

    SAMPLE=$(basename "$BAM" _Aligned.sortedByCoord.out.bam)
    echo "Processing sample: $SAMPLE ..."

    # Detect if paired-end or single-end
    READ_TYPE=$(samtools view -c -f 1 "$BAM")
    if [ "$READ_TYPE" -gt 0 ]; then
      echo "  → Detected paired-end BAM"
      PAIR_FLAG="--paired-end"
    else
      echo "  → Detected single-end BAM"
      PAIR_FLAG=""
    fi

    # Run RSEM quantification
    rsem-calculate-expression \
      --alignments \
      --bam \
      $PAIR_FLAG \
      --estimate-rspd \
      --append-names \
      --output-genome-bam \
      --strandedness reverse \
      --num-threads $THREADS \
      "$BAM" \
      "$REF" \
      "$OUTDIR/$STUDYNAME/${SAMPLE}"

    if [ $? -ne 0 ]; then
      echo "  Error processing $SAMPLE in $STUDYNAME — check log file"
    else
      echo "  Done: $SAMPLE"
    fi
  done
done
