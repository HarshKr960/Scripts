#!/bin/bash

# Set up directories for input and output
DATA_DIR="/home/DATA"
REF_DIR="/home/ref"
OUTPUT_DIR="/home/Output"
QC_DIR="$OUTPUT_DIR/QC"
TRIM_DIR="$OUTPUT_DIR/Trimmed"
ALIGN_DIR="$OUTPUT_DIR/Aligned"
COUNT_DIR="$OUTPUT_DIR/Counts"

# Create output directories if they don't exist
mkdir -p $QC_DIR $TRIM_DIR $ALIGN_DIR $COUNT_DIR

# Reference genome and index
#REFERENCE_GENOME="$REF_DIR/genome.fasta"
HISAT2_INDEX="$REF_DIR/r64/genome"

# Number of threads
THREADS=12

# Start timing the pipeline
START_TIME=$SECONDS

# Step 1: Quality Control (FastQC)
echo "Running FastQC..."
for FILE in $DATA_DIR/*.fastq.gz; do
    fastqc $FILE -t $THREADS -o $QC_DIR
done
echo "FastQC completed."

# Step 2: Trimming (Trimmomatic)
echo "Running Trimmomatic..."
for FILE in $DATA_DIR/*_1.fastq.gz; do
    BASE=$(basename $FILE "_1.fastq.gz")
    TRIMMED_1="$TRIM_DIR/${BASE}_1_trimmed.fastq.gz"
    TRIMMED_2="$TRIM_DIR/${BASE}_2_trimmed.fastq.gz"
    TRIMMOMATIC_ADAPTERS="/root/miniconda3/share/trimmomatic/adapters/NexteraPE-PE.fa"  # Update this path with adapter sequences
    
    trimmomatic PE -threads $THREADS \
   $DATA_DIR/${BASE}_1.fastq.gz $DATA_DIR/${BASE}_2.fastq.gz \
    $TRIMMED_1 /dev/null \
    $TRIMMED_2 /dev/null \
    ILLUMINACLIP:$TRIMMOMATIC_ADAPTERS:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36
done
echo "Trimmomatic completed."


# Step 3: Alignment (HISAT2)
echo "Running HISAT2..."

# Check if trimmed files are present
if [ -z "$(ls -A $TRIM_DIR/*_1_trimmed.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No trimmed files found in $TRIM_DIR"
    exit 1
fi

# Loop over trimmed FASTQ files
for FILE in $TRIM_DIR/*_1_trimmed.fastq.gz; do
    BASE=$(basename $FILE "_1_trimmed.fastq.gz")
    SAM_FILE="$ALIGN_DIR/${BASE}.sam"
    SUMMARY_FILE="$ALIGN_DIR/${BASE}_alignment_summary.txt"

    # Run HISAT2 and save the alignment summary
    echo "Aligning $BASE..."
    if hisat2 -p $THREADS -x $HISAT2_INDEX \
        -1 $TRIM_DIR/${BASE}_1_trimmed.fastq.gz \
        -2 $TRIM_DIR/${BASE}_2_trimmed.fastq.gz \
        -S $SAM_FILE 2> $SUMMARY_FILE; then
        echo "HISAT2 alignment for $BASE completed successfully."
    else
        echo "Error: HISAT2 alignment for $BASE failed. Check $SUMMARY_FILE for details." >&2
        continue  # Skip to the next file if HISAT2 fails
    fi
done

echo "HISAT2 alignment completed."

# Step 4: Convert SAM to BAM and sort (SAMtools)
echo "Converting SAM to sorted BAM..."
for SAM_FILE in $ALIGN_DIR/*.sam; do
    BASE=$(basename $SAM_FILE ".sam")
    BAM_FILE="${ALIGN_DIR}/${BASE}.bam"
    SORTED_BAM_FILE="${ALIGN_DIR}/${BASE}_sorted.bam"

    echo "Processing $SAM_FILE..."

    # Convert SAM to BAM
    if ! samtools view -@ $THREADS -bS $SAM_FILE > $BAM_FILE; then
        echo "Error: Failed to convert $SAM_FILE to BAM" >&2
        continue  # Skip to next file if conversion fails
    fi

    # Sort BAM file
    if ! samtools sort -@ $THREADS $BAM_FILE -o $SORTED_BAM_FILE; then
        echo "Error: Failed to sort $BAM_FILE" >&2
        continue  # Skip to next file if sorting fails
    fi

    # Index sorted BAM file
    if ! samtools index $SORTED_BAM_FILE; then
        echo "Error: Failed to index $SORTED_BAM_FILE" >&2
        continue  # Skip to next file if indexing fails
    fi

    # Cleanup intermediate files
    rm $SAM_FILE $BAM_FILE
done
echo "SAM to BAM conversion and sorting completed."

# Step 5: Feature Counts (featureCounts)
echo "Running featureCounts..."
ANNOTATION_FILE="$REF_DIR/Saccharomyces_cerevisiae.R64-1-1.113.gtf"  # Specify the path to the annotation file (GTF)

featureCounts -T $THREADS -a $ANNOTATION_FILE -o $COUNT_DIR/counts.txt -p $ALIGN_DIR/*_sorted.bam
echo "featureCounts completed."

# Define the input and output file names
INPUT_FILE="$COUNT_DIR/counts.txt"
OUTPUT_FILE="$COUNT_DIR/filtered_counts.txt"

# Notify the user that filtering is starting
echo "Filtering featureCounts output..."

# Run the awk command to filter the count matrix
awk 'BEGIN { OFS="\t" } NR==1 { print; next } { printf "%s\t", $1; for(i=7; i<=NF; i++) printf "%s%s", $i, (i==NF ? "\n" : OFS) }' $INPUT_FILE > $OUTPUT_FILE

# Notify the user when filtering is complete
echo "Filtering complete. Filtered counts saved to $OUTPUT_FILE."

echo "Pipeline completed successfully."
# Calculate total time taken
ELAPSED_TIME=$(( SECONDS - START_TIME ))
echo "Total time taken: $ELAPSED_TIME seconds."