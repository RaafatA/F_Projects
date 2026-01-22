#!/bin/bash

REF="/home/2375894/radwa-scratch/GEO_GTEx/refrence/rsem_reference"
OUTDIR="/home/2375894/radwa-scratch/GEO_GTEx/rsem_results"
THREADS=8
MAX_FRAG_LEN=1000

shopt -s nullglob

for STUDY in /home/2375894/radwa-scratch/GEO_GTEx/mapped_GSE*/; do
  STUDYNAME=$(basename "$STUDY")
  mkdir -p "$OUTDIR/$STUDYNAME"

  echo "========== Processing Study: $STUDYNAME =========="

  for BAM in "$STUDY"/*.bam; do
    [ -e "$BAM" ] || { echo "No BAM files in $STUDYNAME"; continue; }

    SAMPLE=$(basename "$BAM" .Aligned.toTranscriptome.out.bam)
    echo "Processing sample: $SAMPLE ..."

    # Detect paired-end vs single-end
    READ_TYPE=$(samtools view -c -f 1 "$BAM")
    if [ "$READ_TYPE" -gt 0 ]; then
      echo "  → Detected paired-end BAM"
      PAIR_FLAG="--paired-end"
    else
      echo "  → Detected single-end BAM"
      PAIR_FLAG=""
    fi

    # RSEM quantification (parameter-matched)
    rsem-calculate-expression \
      --bam \
      $PAIR_FLAG \
      --estimate-rspd \
      --fragment-length-max $MAX_FRAG_LEN \
      --no-bam-output \
      --num-threads $THREADS \
      "$BAM" \
      "$REF" \
      "$OUTDIR/$STUDYNAME/${SAMPLE}.rsem"

    if [ $? -ne 0 ]; then
      echo "  Error processing $SAMPLE in $STUDYNAME"
    else
      echo "  Done: $SAMPLE"
    fi
  done
done
