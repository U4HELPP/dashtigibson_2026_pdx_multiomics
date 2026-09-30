#!/bin/bash

# Script Name: link_jax_results.sh
# Description: Automates the creation of symbolic links for JAX pipeline results
#              and then runs a modified merge script.
#              It takes a runID (directory name) as an argument, navigates to
#              specific target directories, and creates symbolic links back
#              to the original human and mouse BAM, gene expression, and
#              isoform expression files within the specified runID directory.
#              After linking, it copies, modifies, and executes a data merging script.

# --- Configuration Variables ---
# Base directory where consolidated results links will be created and merge scripts stored.
RESULTS_DIST_BASE="/path/to/bulkRNASeq/jax_pipeline_results"
# Base directory where the original JAX PDX pipeline results are located.
SOURCE_JAXPDX_BASE="/path/to/bulkRNASeq/jaxpdxnf/jaxpdx_results"
# Path to the original merge script template.
MERGE_SCRIPT_SOURCE_PATH="${RESULTS_DIST_BASE}/scripts/merge_rsem_data_250108.sh"

# --- Function to display script usage ---
function display_usage {
    echo "Usage: $0 <runID>"
    echo "  <runID>: The name of the run directory within ${SOURCE_JAXPDX_BASE} "
    echo "           (e.g., 250102_SUS20241011056_results)"
    echo ""
    echo "Example: $0 250102_SUS20241011056_results"
    echo ""
    echo "Note: This script will automatically update the date in the merge script's"
    echo "      RunID (-r argument) and log file names to the current date."
}

# --- Argument Parsing ---
# Check if a runID argument was provided.
if [ -z "$1" ]; then
    echo "Error: runID argument is missing."
    display_usage
    exit 1
fi

# Assign the first command-line argument to RUN_ID.
RUN_ID="$1"
# Construct the full path to the source results directory for this specific run.
FULL_SOURCE_RUN_PATH="${SOURCE_JAXPDX_BASE}/${RUN_ID}"

# --- Validate Source Directory ---
# Check if the full source directory exists.
if [ ! -d "${FULL_SOURCE_RUN_PATH}" ]; then
    echo "Error: Source directory '${FULL_SOURCE_RUN_PATH}' does not exist."
    echo "Please ensure the runID is correct and the directory exists."
    exit 1
fi

echo "--- Starting JAX Pipeline Automation for RunID: ${RUN_ID} ---"
echo "Source path for linking: ${FULL_SOURCE_RUN_PATH}"
echo "Target base directory: ${RESULTS_DIST_BASE}"
echo ""

# --- Create Symbolic Links ---
# Loop through human and mouse species.
for species in "human" "mouse"; do
    # Determine the file pattern based on species for gene and isoform expression.
    FILE_PATTERN_GENE="*${species}.genes.results"
    FILE_PATTERN_ISOFORM="*${species}.isoforms.results"

    echo "Processing ${species} files..."

    # --- BAM Files ---
    TARGET_DIR="${RESULTS_DIST_BASE}/${species}_bam"
    SOURCE_PATH_BAM="${FULL_SOURCE_RUN_PATH}/*_${species}/bam/*"
    echo "  Navigating to: ${TARGET_DIR}"
    if [ ! -d "${TARGET_DIR}" ]; then
        echo "  Warning: Target directory '${TARGET_DIR}' does not exist. Creating it."
        mkdir -p "${TARGET_DIR}" || { echo "Error: Failed to create directory ${TARGET_DIR}"; exit 1; }
    fi
    # Use pushd/popd to safely navigate and return to the original directory.
    pushd "${TARGET_DIR}" > /dev/null
    echo "  Linking BAM files from: ${SOURCE_PATH_BAM}"
    # The 'shopt -s nullglob' ensures that if no files match the pattern,
    # the wildcard expands to nothing, preventing 'ln: cannot access' errors.
    shopt -s nullglob
    ln -sf ${SOURCE_PATH_BAM} .
    shopt -u nullglob # Turn off nullglob
    popd > /dev/null
    echo "  BAM linking complete for ${species}."

    # --- Gene Expression Files ---
    TARGET_DIR="${RESULTS_DIST_BASE}/${species}_gene_exp"
    SOURCE_PATH_GENE_EXP="${FULL_SOURCE_RUN_PATH}/*_${species}/${FILE_PATTERN_GENE}"
    echo "  Navigating to: ${TARGET_DIR}"
    if [ ! -d "${TARGET_DIR}" ]; then
        echo "  Warning: Target directory '${TARGET_DIR}' does not exist. Creating it."
        mkdir -p "${TARGET_DIR}" || { echo "Error: Failed to create directory ${TARGET_DIR}"; exit 1; }
    fi
    pushd "${TARGET_DIR}" > /dev/null
    echo "  Linking Gene Expression files from: ${SOURCE_PATH_GENE_EXP}"
    shopt -s nullglob
    ln -sf ${SOURCE_PATH_GENE_EXP} .
    shopt -u nullglob
    popd > /dev/null
    echo "  Gene Expression linking complete for ${species}."

    # --- Isoform Expression Files ---
    TARGET_DIR="${RESULTS_DIST_BASE}/${species}_isoform_exp"
    SOURCE_PATH_ISOFORM_EXP="${FULL_SOURCE_RUN_PATH}/*_${species}/${FILE_PATTERN_ISOFORM}"
    echo "  Navigating to: ${TARGET_DIR}"
    if [ ! -d "${TARGET_DIR}" ]; then
        echo "  Warning: Target directory '${TARGET_DIR}' does not exist. Creating it."
        mkdir -p "${TARGET_DIR}" || { echo "Error: Failed to create directory ${TARGET_DIR}"; exit 1; }
    fi
    pushd "${TARGET_DIR}" > /dev/null
    echo "  Linking Isoform Expression files from: ${SOURCE_PATH_ISOFORM_EXP}"
    shopt -s nullglob
    ln -sf ${SOURCE_PATH_ISOFORM_EXP} .
    shopt -u nullglob
    popd > /dev/null
    echo "  Isoform Expression linking complete for ${species}."
    echo ""
done

echo "--- All symbolic links created successfully for RunID: ${RUN_ID} ---"
echo ""

# --- Prepare and Run Merge Script ---

# Get current date in YY.MM.DD format for -r argument and log files.
CURRENT_DATE_DOTTED=$(date +%y.%m.%d)
# Get current date in YYMMDD format for the new merge script filename.
CURRENT_DATE_NOSLASH=$(date +%y%m%d)

# Define the new filename for the copied merge script.
MERGE_SCRIPT_TARGET_FILENAME="merge_rsem_data_${CURRENT_DATE_NOSLASH}.sh"
MERGE_SCRIPT_TARGET_PATH="${RESULTS_DIST_BASE}/scripts/${MERGE_SCRIPT_TARGET_FILENAME}"

echo "--- Preparing and Running Merge Script ---"

# 1. Check if the source merge script exists.
if [ ! -f "${MERGE_SCRIPT_SOURCE_PATH}" ]; then
    echo "Error: Source merge script '${MERGE_SCRIPT_SOURCE_PATH}' not found."
    echo "Please ensure the original merge script exists at this location."
    exit 1
fi

# 2. Copy the merge script to the dated filename.
echo "Copying merge script from '${MERGE_SCRIPT_SOURCE_PATH}' to '${MERGE_SCRIPT_TARGET_PATH}'"
cp "${MERGE_SCRIPT_SOURCE_PATH}" "${MERGE_SCRIPT_TARGET_PATH}" || { echo "Error: Failed to copy merge script."; exit 1; }

# 3. Edit the copied merge script to update RunID and log file dates.
echo "Updating RunID (-r argument) and log file names in '${MERGE_SCRIPT_TARGET_FILENAME}'"

# Replace the hardcoded date '25.01.08_' in the -r argument with the current date.
# This ensures the new merged files will have the correct date prefix.
sed -i "s/25\.01\.08_/${CURRENT_DATE_DOTTED}_/g" "${MERGE_SCRIPT_TARGET_PATH}"

# Replace the hardcoded date '25.01.08_human' and '25.01.08_mouse' in log file names with the current date.
sed -i "s/25\.01\.08_human/${CURRENT_DATE_DOTTED}_human/g" "${MERGE_SCRIPT_TARGET_PATH}"
sed -i "s/25\.01\.08_mouse/${CURRENT_DATE_DOTTED}_mouse/g" "${MERGE_SCRIPT_TARGET_PATH}"

# 4. Make the copied script executable.
chmod +x "${MERGE_SCRIPT_TARGET_PATH}" || { echo "Error: Failed to make merge script executable."; exit 1; }

# 5. Run the merge script.
echo "Executing the generated merge script: '${MERGE_SCRIPT_TARGET_FILENAME}'"
# Use 'source' as requested, to ensure environment variables (like module loads) are preserved.
source "${MERGE_SCRIPT_TARGET_PATH}"

echo "--- Merge script execution complete ---"
echo "--- JAX Pipeline Automation finished ---"