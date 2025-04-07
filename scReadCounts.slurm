#!/bin/bash
#SBATCH --job-name=scReadCounts        # Job name
#SBATCH --output=logs/scReadCounts.%A_%a.out # Standard output log (%A=jobID, %a=array task ID)
#SBATCH --error=logs/scReadCounts.%A_%a.err  # Standard error log
#SBATCH --nodes=1                      # Run all tasks on a single node
#SBATCH --ntasks=1                     # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem-per-cpu=4gb              # Memory per CPU core
#SBATCH --time=20:00:00                # Time limit hrs:min:sec
#SBATCH --qos=your_qos                 # Specify the QoS partition (e.g., u54care2)
#SBATCH --account=your_account         # Specify the account (e.g., u54care2)
#SBATCH --array=1-N                    # Set N to the total number of samples in your run_list.txt

### Purpose: Run SCReadCounts on single-cell BAM files using pre-split mitochondrial SNV lists.
### Input:   - Cell Ranger output BAM file (possorted_genome_bam.bam)
###          - Cell Ranger filtered barcodes file (barcodes.tsv)
###          - A run list file specifying project directories and run names
###          - Directory containing pre-split mitochondrial SNV lists
### Output:  - SCReadCounts output directory containing multiple TSV files (one per SNV chunk).

# --- Configuration ---
RUN_LIST_FILE="/path/to/your/run_list.txt"             # IMPORTANT: Update this path
LOG_DIR="logs"                                         # Directory for SLURM output/error logs
SCREADCOUNTS_BIN_DIR="/path/to/SCReadCounts/bin/dir"   # IMPORTANT: Update path to SCReadCounts binary directory
SNV_LIST_DIR="/path/to/your/split_snv_lists/400bases"  # IMPORTANT: Update path to dir with mt_vars_range_M_*.txt files
INTERVAL=400                                           # Interval size used for splitting SNV lists
NUM_INTERVALS=42                                       # Number of split SNV files (0 to 41)

# --- Create Log Directory ---
mkdir -p $LOG_DIR

# --- Environment Setup ---
# SCReadCounts might be a standalone binary, no conda env needed unless specified
echo "Setting up environment..."
# Ensure the binary directory is executable or add it to PATH if needed
# export PATH=$SCREADCOUNTS_BIN_DIR:$PATH

# Check if scReadCounts command is available/executable
SCREADCOUNTS_EXE="${SCREADCOUNTS_BIN_DIR}/scReadCounts" # Adjust executable name if different
if [ ! -x "$SCREADCOUNTS_EXE" ]; then
    echo "Error: SCReadCounts executable not found or not executable: $SCREADCOUNTS_EXE"
    exit 1
fi
echo "Using SCReadCounts executable: ${SCREADCOUNTS_EXE}"

# --- Get Sample Info from Run List ---
TASK_ID=$SLURM_ARRAY_TASK_ID
echo "Processing SLURM Array Task ID: $TASK_ID"

# Check if RUN_LIST_FILE exists
if [ ! -f "$RUN_LIST_FILE" ]; then
    echo "Error: Run list file not found: $RUN_LIST_FILE"
    exit 1
fi

# Read sample information (adjust column numbers if needed)
PROJECT_DIR=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $2}' "$RUN_LIST_FILE")
RUN_NAME=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $3}' "$RUN_LIST_FILE")
# SAMPLE_NAME=$(awk -F'\t' -v id="$TASK_ID" '$1 == id {print $4}' "$RUN_LIST_FILE") # Use if needed

# Validate extracted info
if [ -z "$PROJECT_DIR" ] || [ -z "$RUN_NAME" ]; then
  echo "Error: Could not parse required information (ProjectDir, RunName) for Task ID $TASK_ID from $RUN_LIST_FILE"
  exit 1
fi

# --- Define Input/Output Paths ---
CELLRANGER_OUT_DIR="${PROJECT_DIR}/cellranger_RNA_out/${RUN_NAME}_RNA/outs"
BAM_FILE="${CELLRANGER_OUT_DIR}/possorted_genome_bam.bam"
BARCODE_DIR="${CELLRANGER_OUT_DIR}/filtered_feature_bc_matrix"
BARCODE_FILE_GZ="${BARCODE_DIR}/barcodes.tsv.gz"
BARCODE_FILE="${BARCODE_DIR}/barcodes.tsv" # Unzipped version path
OUTPUT_DIR="${PROJECT_DIR}/scReadCounts_out/${RUN_NAME}/MT" # Output directory for this sample
OUTPUT_FILE_PREFIX="${OUTPUT_DIR}/scReadCounts"            # Base name for output chunk files

# --- Input File Checks ---
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi
if [ ! -f "$BARCODE_FILE_GZ" ]; then
    echo "Error: Gzipped barcode file not found: $BARCODE_FILE_GZ"
    exit 1
fi
if [ ! -d "$SNV_LIST_DIR" ]; then
    echo "Error: SNV list directory not found: $SNV_LIST_DIR"
    exit 1
fi

# --- Prepare Barcodes ---
# SCReadCounts requires an unzipped barcode file.
# Warning: Modifies Cell Ranger output directory. Consider copying first.
echo "Unzipping barcode file: ${BARCODE_FILE_GZ}"
gunzip -f "${BARCODE_FILE_GZ}"
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
echo "Starting SCReadCounts Analysis"
echo "SLURM Job ID:      $SLURM_JOB_ID"
echo "SLURM Task ID:     $TASK_ID"
echo "Run Name:          $RUN_NAME"
echo "Project Dir:       $PROJECT_DIR"
echo "Input BAM:         $BAM_FILE"
echo "Input Barcodes:    $BARCODE_FILE"
echo "SNV List Dir:      $SNV_LIST_DIR"
echo "Output Dir:        $OUTPUT_DIR"
echo "CPUs:              $SLURM_CPUS_PER_TASK"
echo "----------------------------------------"

# SCReadCounts expects to be run from the directory containing the barcode file
echo "Changing directory to barcode directory: ${BARCODE_DIR}"
cd "$BARCODE_DIR" || { echo "Error: Failed to change directory to $BARCODE_DIR"; exit 1; }

echo "Running SCReadCounts in chunks..."
overall_status=0
for i in $(seq 0 $((NUM_INTERVALS - 1)))
do
    start=$((i * INTERVAL))
    end=$((start + INTERVAL - 1))
    mt_vars_chunk="${SNV_LIST_DIR}/mt_vars_range_M_${start}_${end}.txt"
    output_chunk="${OUTPUT_FILE_PREFIX}_${start}.tsv"

    if [ ! -f "$mt_vars_chunk" ]; then
        echo "Warning: SNV chunk file not found, skipping: $mt_vars_chunk"
        continue
    fi

    echo "Processing chunk: ${mt_vars_chunk} -> ${output_chunk}"

    # Run SCReadCounts for the chunk
    "${SCREADCOUNTS_EXE}" -r "${BAM_FILE}" \
                         -s "${mt_vars_chunk}" \
                         -C CellRanger \
                         -b "${BARCODE_FILE}" \
                         --umicount CellRanger \
                         --directional True \
                         -t "${SLURM_CPUS_PER_TASK}" \
                         --minreads=10 \
                         -o "${output_chunk}"

    # Check exit status for the chunk
    status=$?
    if [ $status -ne 0 ]; then
      echo "Error: SCReadCounts failed for chunk ${start}-${end} with exit status $status"
      overall_status=1 # Mark that at least one chunk failed
    fi

    # Sleep briefly - purpose unclear (filesystem load?), keep as in original script
    echo "Sleeping for 60 seconds..."
    sleep 60
done

# Go back to the original submission directory (optional)
# cd -

# Final status check
if [ $overall_status -ne 0 ]; then
  echo "Error: One or more SCReadCounts chunks failed for $RUN_NAME (Task ID $TASK_ID)."
  # Optional: Re-zip barcode file
  # gzip "${BARCODE_FILE}"
  exit 1
else
  echo "SCReadCounts chunk processing completed for $RUN_NAME (Task ID $TASK_ID)."
  echo "Output chunks generated in: $OUTPUT_DIR"
fi

# Optional: Re-zip the barcode file
# echo "Re-zipping barcode file: ${BARCODE_FILE}"
# gzip "${BARCODE_FILE}"

echo "----------------------------------------"
echo "SLURM task $TASK_ID finished."

exit 0