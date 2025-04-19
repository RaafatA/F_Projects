#!/bin/bash
# Exit immediately if a command exits with a non-zero status.
set -e

#############################
# 1. Create and activate a new conda environment
#############################
# Create a new environment (here named "miRDeep2_env")
conda create -n miRDeep2_env -y
conda activate miRDeep2_env

# Install required packages from Bioconda:
# - miRDeep2 (includes mapper.pl and miRDeep2.pl)
# - sra-tools (for downloading SRA data)
# - fastqc (for quality control)
# - cutadapt (for adapter trimming)
# - bowtie (for read alignment)
# - samtools (for additional processing)
conda install -c bioconda mirdeep2 sra-tools fastqc cutadapt bowtie samtools -y

#############################
# 2. Download sequencing data from SRA
#############################
# Replace "SRR1234567" with your actual SRA accession number.
# Here we download and split paired-end (or single-end) data.
fastq-dump --split-files SRR1234567

#############################
# 3. Pre-process raw reads
#############################
# Run FastQC for an initial quality check
fastqc SRR1234567_1.fastq

# Trim adapter sequences with cutadapt.
# (The adapter sequence below is a common miRNA-seq adapter; adjust if needed.)
cutadapt -a TGGAATTCTCGGGTGCCAAGG -o SRR1234567_trimmed.fastq SRR1234567_1.fastq

# Optionally, re-run FastQC on the trimmed file to verify quality improvement
fastqc SRR1234567_trimmed.fastq

#############################
# 4. Download and index the reference genome
#############################
# For novel miRNA discovery, mapping is done against the full human reference.
# Download the human reference genome (here using Ensembl’s GRCh38 primary assembly)
wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Build the Bowtie index for the human genome (required by miRDeep2 mapper)
bowtie-build Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38_index

#############################
# 5. Download miRBase reference sequences for known miRNAs
#############################
# For quantifying known miRNAs, you need the mature and hairpin sequences.
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
gunzip mature.fa.gz
gunzip hairpin.fa.gz

# Extract only the human (hsa) sequences from the downloaded files.
# (Assumes headers start with ">hsa-")
grep -A1 '^>hsa-' mature.fa > hsa.mature.fa
grep -A1 '^>hsa-' hairpin.fa > hsa.hairpin.fa

#############################
# 6. Prepare reads for miRDeep2 analysis (using mapper.pl)
#############################
# mapper.pl collapses identical reads and maps them to the genome.
# Parameters used:
#   -e : remove reads with non-canonical letters
#   -h : remove low complexity reads
#   -m : collapse reads into a fasta file
#   -j : remove reads with counts of 1 if desired (optional)
#   -l 18 : discard reads shorter than 18 nucleotides
#   -s : output file for collapsed reads
#   -t : output mapping results in .arf format (needed by miRDeep2.pl)
#   -p : prefix of the Bowtie index built earlier
mapper.pl SRR1234567_trimmed.fastq -e -d -h -m -j -l 18 -s reads_collapsed.fa -t reads_vs_genome.arf -p hg38_index

#############################
# 7. Run miRDeep2 for miRNA quantification and discovery
#############################
# miRDeep2.pl takes the collapsed reads, the genome, mapping file, and the miRBase references.
# Note: the second parameter for mature sequences is repeated – this allows using the same file for both known and candidate quantification.
awk 'print $1' Homo_sapiens.GRCh38.dna.primary_assembly.fa > one_column_ref.fa
miRDeep2.pl reads_collapsed.fa one_column_ref.fa reads_vs_genome.arf hsa.mature.fa hsa.hairpin.fa -t hsa

# The output directory will contain files including:
#   - HTML reports with predicted miRNAs and their scores
#   - Expression quantification files (often a count matrix for known miRNAs)
# You can inspect these outputs (for example, the "expression_analyses" folder) to retrieve your count matrix.
