#!/bin/bash
STAR --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir /home/2375894/radwa-scratch/GEO_GTEx/refrence \
    --genomeFastaFiles /home/2375894/radwa-scratch/GEO_GTEx/refrence/GRCh38.primary_assembly.genome.fa.gz \
    --sjdbGTFfile /home/2375894/radwa-scratch/GEO_GTEx/refrence/gencode.v49.annotation.gtf.gz \
    --limitGenomeGenerateRAM 64000000000 \
