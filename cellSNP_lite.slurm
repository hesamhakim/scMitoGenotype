#!/bin/bash
#SBATCH --job-name=cellSNP_lite      # Job name
#SBATCH --output=logs/cellSNP_lite.%A_%a.out # Standard output log (%A=jobID, %a=array task ID)
#SBATCH --error=logs/cellSNP_lite.%A_%a.err  # Standard error log
#SBATCH --nodes=1                    # Run all tasks on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=16           # Number of CPU cores per task
#SBATCH --mem-per-cpu=8gb            # Memory per CPU core
#SBATCH --time=30:00:00              # Time limit hrs:min:sec
#SBATCH --qos=your_qos               # Specify the QoS partition (e.g., u54care2)
#SBATCH --account=your_account       # Specify the account (e.g., u54care2)
#SBATCH --array=1-N                  # Set N to the total number of samples in your run_list.txt

### Purpose: Run cellsnp-lite on single-cell BAM files to genotype mitochondrial variants.
### Input:   - Cell Ranger output BAM file (possorted_genome_bam.bam)
###          - Cell Ranger filtered barcodes file (barcodes.tsv)
###          - A run list file specifying project directories and run names
### Output:  - cellsnp-lite output directory containing VCF and MTX files.

# --- Configuration ---
RUN_LIST_FILE="/path/to/your/run_list.txt" # IMPORTANT: Update this path
LOG_DIR="logs"                             # Directory for SLURM output/error logs
CONDA_ENV_NAME="mquad"                     # Name of the conda environment with cellsnp-lite

# --- Create Log Directory ---
mkdir -p $LOG_DIR

# --- Environment Setup ---
echo "Setting up environment..."
source ~/.bashrc # Source your bash profile if needed
conda activate "${CONDA_ENV_NAME}"
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment '${CONDA_ENV_NAME}'"
    exit 1
fi

# Check if cellsnp-lite command is available
if ! command -v cellsnp-lite &> /dev/null; then
    echo "Error: 'cellsnp-lite' command not found in environment '${CONDA_ENV_NAME}'."
    exit 1
fi
echo "Using cellsnp-lite from environment: ${CONDA_ENV_NAME}"

# --- Get Sample Info from Run List ---
TASK_ID=$SLURM_ARRAY_TASK_ID
echo "Processing SLURM Array Task ID: $TASK_ID"

# Check if RUN_LIST_FILE exists
if [ ! -f "$RUN_LIST_FILE" ]; then
    echo "Error: Run list file not found: $RUN_LIST_FILE"
    exit 1
fi

# Read sample information from the run list based on the task ID
# Assumes run_list.txt is TAB-delimited with columns:
# 1: ID (numeric, matching SLURM array index)
# 2: ProjectDir (base directory for input/output)
# 3: RunName (used in Cell Ranger output path)
# 4: SampleName (optional, potentially same as RunName)
# --- Adjust column numbers ($2, $3, etc.) if your run_list.txt format differs ---
PROJECT_DIR=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $2}' "$RUN_LIST_FILE")
RUN_NAME=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $3}' "$RUN_LIST_FILE")
# SAMPLE_NAME=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $4}' "$RUN_LIST_FILE") # Use if needed

# Validate extracted info
if [ -z "$PROJECT_DIR" ] || [ -z "$RUN_NAME" ]; then
  echo "Error: Could not parse required information (ProjectDir, RunName) for Task ID $TASK_ID from $RUN_LIST_FILE"
  echo "Check line $TASK_ID in the file."
  exit 1
fi

# --- Define Input/Output Paths ---
CELLRANGER_OUT_DIR="${PROJECT_DIR}/cellranger_RNA_out/${RUN_NAME}_RNA/outs"
BAM_FILE="${CELLRANGER_OUT_DIR}/possorted_genome_bam.bam"
BARCODE_DIR="${CELLRANGER_OUT_DIR}/filtered_feature_bc_matrix"
BARCODE_FILE_GZ="${BARCODE_DIR}/barcodes.tsv.gz"
BARCODE_FILE="${BARCODE_DIR}/barcodes.tsv" # Unzipped version path
OUTPUT_DIR="${PROJECT_DIR}/scSNP_lite_out/${RUN_NAME}"

# --- Input File Checks ---
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi
if [ ! -f "$BARCODE_FILE_GZ" ]; then
    echo "Error: Gzipped barcode file not found: $BARCODE_FILE_GZ"
    exit 1
fi

# --- Prepare Barcodes ---
# cellsnp-lite requires an unzipped barcode file.
# Warning: This modifies the Cell Ranger output directory.
# Consider copying the barcode file elsewhere before unzipping if modification is undesirable.
echo "Unzipping barcode file: ${BARCODE_FILE_GZ}"
gunzip -f "${BARCODE_FILE_GZ}" # Use -f to force overwrite if .tsv already exists
if [ $? -ne 0 ]; then
    echo "Error: Failed to unzip barcode file: ${BARCODE_FILE_GZ}"
    exit 1
fi
if [ ! -f "$BARCODE_FILE" ]; then
    echo "Error: Unzipped barcode file not found after gunzip attempt: $BARCODE_FILE"
    exit 1
fi

# --- Create Output Directory ---
mkdir -p "${OUTPUT_DIR}"
if [ $? -ne 0 ]; then
    echo "Error: Failed to create output directory: ${OUTPUT_DIR}"
    exit 1
fi

# --- Execution ---
echo "----------------------------------------"
echo "Starting cellsnp-lite Analysis"
echo "SLURM Job ID:      $SLURM_JOB_ID"
echo "SLURM Task ID:     $TASK_ID"
echo "Run Name:          $RUN_NAME"
echo "Project Dir:       $PROJECT_DIR"
echo "Input BAM:         $BAM_FILE"
echo "Input Barcodes:    $BARCODE_FILE"
echo "Output Directory:  $OUTPUT_DIR"
echo "CPUs:              $SLURM_CPUS_PER_TASK"
echo "----------------------------------------"

# Run cellsnp-lite for mitochondrial genome (chrM)
# Mode 1a: droplet based single cells without given SNPs (generates VCF)
cellsnp-lite -s "$BAM_FILE" \
             -b "$BARCODE_FILE" \
             -O "$OUTPUT_DIR" \
             --chrom M \
             --UMItag UB \
             --minMAF 0.01 \
             --minCOUNT 10 \
             --genotype \
             --gzip \
             -p "$SLURM_CPUS_PER_TASK"

# Check exit status
status=$?
if [ $status -ne 0 ]; then
  echo "Error: cellsnp-lite failed for $RUN_NAME (Task ID $TASK_ID) with exit status $status"
  # Optional: Re-zip the barcode file if needed by other processes
  # gzip "${BARCODE_FILE}"
  exit $status
else
  echo "cellsnp-lite completed successfully for $RUN_NAME (Task ID $TASK_ID)."
  echo "Output generated in: $OUTPUT_DIR"
fi

# Optional: Re-zip the barcode file if needed by other processes
# echo "Re-zipping barcode file: ${BARCODE_FILE}"
# gzip "${BARCODE_FILE}"

echo "----------------------------------------"
echo "SLURM task $TASK_ID finished."

exit 0