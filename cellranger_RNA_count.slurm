#!/bin/bash
#SBATCH --job-name=cellranger_count  # Job name
#SBATCH --output=logs/cellranger_%A_%a.out # Standard output log (%A=jobID, %a=array task ID)
#SBATCH --error=logs/cellranger_%A_%a.err  # Standard error log
#SBATCH --nodes=1                    # Run all tasks on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=16           # Number of CPU cores per task (adjust as needed)
#SBATCH --mem=64G                    # Total memory for the job (adjust as needed, e.g., 128G)
#SBATCH --time=24:00:00              # Time limit hrs:min:sec (adjust as needed)
#SBATCH --qos=your_qos               # Specify the QoS partition, e.g., u54care2
#SBATCH --account=your_account       # Specify the account, e.g., u54care2
#SBATCH --array=1-N                  # Set N to the total number of samples in your run_list.txt

### Purpose: Run 10x Cell Ranger 'count' pipeline on FASTQ files using SLURM arrays.
### Input:   - FASTQ files for multiple samples
###          - Cell Ranger reference transcriptome
###          - A run list file specifying sample details
### Output:  - Cell Ranger output directory for each sample

# --- Configuration ---
RUN_LIST_FILE="/path/to/your/run_list.txt" # IMPORTANT: Update this path
LOG_DIR="logs"                             # Directory for SLURM output/error logs

# --- Create Log Directory ---
mkdir -p $LOG_DIR

# --- Environment Setup ---
echo "Setting up environment..."
# Example: Load Cell Ranger module (adjust version as needed)
# module load cellranger/7.1.0

# Example: Activate Conda environment (if you installed Cell Ranger via Conda)
# source ~/.bashrc
# conda activate cellranger_env

# Check if cellranger command is available
if ! command -v cellranger &> /dev/null; then
    echo "Error: 'cellranger' command not found. Please ensure it's in your PATH or load the necessary module/environment."
    exit 1
fi
echo "Using Cell Ranger version: $(cellranger --version)"

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
# 2: SampleName (prefix for FASTQs, used for --sample)
# 3: FastqDir (path to directory containing FASTQ files for this sample)
# 4: TranscriptomePath (path to Cell Ranger reference)
# 5: ExpectedCells (number, for --expect-cells)
# 6: ProjectDir (base directory for output)
# --- Adjust column numbers ($2, $3, etc.) if your run_list.txt format differs ---
SAMPLE_NAME=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $2}' "$RUN_LIST_FILE")
FASTQ_DIR=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $3}' "$RUN_LIST_FILE")
TRANSCRIPTOME_PATH=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $4}' "$RUN_LIST_FILE")
EXPECTED_CELLS=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $5}' "$RUN_LIST_FILE")
PROJECT_DIR=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $6}' "$RUN_LIST_FILE")

# Validate extracted info
if [ -z "$SAMPLE_NAME" ] || [ -z "$FASTQ_DIR" ] || [ -z "$TRANSCRIPTOME_PATH" ] || [ -z "$EXPECTED_CELLS" ] || [ -z "$PROJECT_DIR" ]; then
  echo "Error: Could not parse complete information for Task ID $TASK_ID from $RUN_LIST_FILE"
  echo "Check line $TASK_ID in the file and ensure all columns are present."
  exit 1
fi
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: FASTQ directory not found for Task ID $TASK_ID: $FASTQ_DIR"
    exit 1
fi
if [ ! -d "$TRANSCRIPTOME_PATH" ]; then
    echo "Error: Transcriptome directory not found for Task ID $TASK_ID: $TRANSCRIPTOME_PATH"
    exit 1
fi

# --- Define Output and Run ID ---
# Create a unique ID for this specific Cell Ranger run
RUN_ID="${SAMPLE_NAME}_${SLURM_ARRAY_JOB_ID}_${TASK_ID}"
# Define the parent directory where Cell Ranger output folders will be created
OUTPUT_PARENT_DIR="${PROJECT_DIR}/cellranger_RNA_out"

# Create the parent output directory if it doesn't exist
mkdir -p "$OUTPUT_PARENT_DIR"

# --- Resource Allocation ---
# Use SLURM variables to inform Cell Ranger
CORES=$SLURM_CPUS_PER_TASK
# Extract numeric value from SLURM memory request (e.g., "64G" -> 64)
# Use SLURM_MEM_PER_NODE as we requested memory for the whole job (--mem=64G)
MEM_GB=$(echo $SLURM_MEM_PER_NODE | grep -o '[0-9]*')

# If you used --mem-per-cpu instead, calculate total memory:
# MEM_PER_CPU_MB=$(echo $SLURM_MEM_PER_CPU | grep -o '[0-9]*') # Assumes MB if no unit like 'G'
# MEM_GB=$(( SLURM_CPUS_PER_TASK * MEM_PER_CPU_MB / 1024 ))

# Handle case where memory extraction might fail
if [ -z "$MEM_GB" ]; then
    echo "Warning: Could not automatically determine memory allocation from SLURM. Defaulting to 8GB for Cell Ranger."
    MEM_GB=8 # Set a fallback default
fi


# --- Execution ---
echo "----------------------------------------"
echo "Starting Cell Ranger Count Analysis"
echo "SLURM Job ID:      $SLURM_JOB_ID"
echo "SLURM Task ID:     $TASK_ID"
echo "Run List Line:     $TASK_ID"
echo "Run ID (--id):     $RUN_ID"
echo "Sample (--sample): $SAMPLE_NAME"
echo "FASTQ Dir:         $FASTQ_DIR"
echo "Transcriptome:     $TRANSCRIPTOME_PATH"
echo "Expected Cells:    $EXPECTED_CELLS"
echo "Project Dir:       $PROJECT_DIR"
echo "Output Parent Dir: $OUTPUT_PARENT_DIR"
echo "Cores (--localcores): $CORES"
echo "Memory (--localmem):  $MEM_GB GB"
echo "----------------------------------------"

# Change to the desired parent output directory before running Cell Ranger
# Cell Ranger will create its output folder (named after --id) inside this directory
cd "$OUTPUT_PARENT_DIR" || { echo "Error: Failed to change directory to $OUTPUT_PARENT_DIR"; exit 1; }

echo "Current working directory: $(pwd)"
echo "Running Cell Ranger..."

# Execute cellranger count
cellranger count --id="${RUN_ID}" \
                 --sample="${SAMPLE_NAME}" \
                 --fastqs="${FASTQ_DIR}" \
                 --transcriptome="${TRANSCRIPTOME_PATH}" \
                 --expect-cells="${EXPECTED_CELLS}" \
                 --localcores="${CORES}" \
                 --localmem="${MEM_GB}"

# Check exit status
status=$?
if [ $status -ne 0 ]; then
  echo "Error: Cell Ranger count failed for $RUN_ID (Task ID $TASK_ID) with exit status $status"
  # cd - # Go back to original directory if needed
  exit $status
else
  echo "Cell Ranger count completed successfully for $RUN_ID (Task ID $TASK_ID)."
  echo "Output generated in: $(pwd)/${RUN_ID}"
fi

# cd - # Go back to original directory if needed

echo "----------------------------------------"
echo "SLURM task $TASK_ID finished."

exit 0