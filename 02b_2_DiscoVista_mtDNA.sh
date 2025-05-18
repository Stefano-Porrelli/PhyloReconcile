#!/bin/bash
#------------------------------------------------------------------------------#
#                           PhyloReconcile                                     #
#                                                                              #
# Script_02b - 2. Visualization of gene tree and species tree discordance      #
# with DiscoVista                                                              #
#------------------------------------------------------------------------------#
BASE_DIR="$(pwd)/PhyloReconcile"
MTDNA_DIR="${BASE_DIR}/02_phylogenies/02a_mtDNA_phylogeny"
MTDNA_FILE="MtDNA_loci.treefile"
GENE_TREES_ML="${MTDNA_DIR}/${MTDNA_FILE}"
CONCAT_TREE="${BASE_DIR}/02_phylogenies/02a_mtDNA_phylogeny/trees/MtDNA_concat.tree"
DISCOVISTA_DIR="${BASE_DIR}/08_DiscoVista_mtDNA"
PARAMETERS="${DISCOVISTA_DIR}/parameters"
DISCOVISTA_TREES="${DISCOVISTA_DIR}/species"
DISCOVISTA_GENES="${DISCOVISTA_DIR}/genetrees"
DISCOVISTA_RESULTS="${DISCOVISTA_DIR}/results"
DISCOVISTA_FREQ="${DISCOVISTA_DIR}/relativeFreq"
DISCOVISTA_INPUTS_CLADEMONOPHYLY="${BASE_DIR}/01_initial_data/Input_files/DiscoVista_parameters_mtDNA/annotation_mtDNA_Clade_Monophyly.txt"
DISCOVISTA_INPUTS_CLADEMONOPHYLY2="${BASE_DIR}/01_initial_data/Input_files/DiscoVista_parameters_mtDNA/annotation_mtDNA_Clade_Monophyly2.txt"

# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate PhyloReconcile

# Create necessary directories
mkdir -p "${DISCOVISTA_DIR}" "${PARAMETERS}" "${DISCOVISTA_TREES}" "${DISCOVISTA_RESULTS}"
# Create directory for individual tree files
mkdir -p "${DISCOVISTA_GENES}"
# Process the single MtDNA_loci.treefile file containing multiple trees
counter=1
while IFS= read -r line
do
    # Skip empty lines
    if [[ -n "$line" ]]; then
        # Create the parent directory structure
        mkdir -p "${DISCOVISTA_GENES}/${counter}/${counter}-ML"
        # Write the tree line to a new file
        echo "$line" > "${DISCOVISTA_GENES}/${counter}/${counter}-ML/estimated_gene_trees.tree"
        echo "Created tree file ${counter}"
        ((counter++))
    fi
done < "${GENE_TREES_ML}"
echo "All trees have been extracted and organized into numbered directories."

# Discordance analysis on gene trees 
# Show the proportion of gene trees for which clades are supported/rejected
cd "${DISCOVISTA_DIR}" || exit
mkdir ./results/Clade_monophyly
# Generate clade-defs file for Tribe/subtribe level (from annotation file)
echo "Generating clade definitions..."
cp "${DISCOVISTA_INPUTS_CLADEMONOPHYLY}" "${DISCOVISTA_DIR}/parameters"
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista generate_clade-defs.py ./parameters/annotation_mtDNA_Clade_Monophyly.txt ./parameters/clade-defs_mtDNA_Clade_Monophyly.txt
# Run DiscoVista analysis
echo "Running DiscoVista analysis..."
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py \
    -c parameters/clade-defs_mtDNA_Clade_Monophyly.txt \
    -p genetrees/ \
    -t 95 \
    -m 1 \
    -o results/Clade_monophyly
echo "DiscoVista analysis completed. Results are in ${DISCOVISTA_RESULTS}"

# Relative frequencies analysis 
# Copy and organize species and gene trees for Relative frequency analysis
cd "${DISCOVISTA_DIR}" || exit
mkdir -p "${DISCOVISTA_FREQ}"/MtDNA_genes-ML
cp "${CONCAT_TREE}" ./species/estimated_species_tree.tree
# retrieve mtDNA species tree and ML gene trees
# ML gene trees need to be on a single file. one tree per line
cp "${DISCOVISTA_TREES}"/estimated_species_tree.tree "${DISCOVISTA_FREQ}"/MtDNA_genes-ML
cp "${GENE_TREES_ML}" "${DISCOVISTA_FREQ}"/MtDNA_genes-ML/estimated_gene_trees.tree
cp "${DISCOVISTA_INPUTS_CLADEMONOPHYLY2}" "${DISCOVISTA_DIR}/parameters"

# Relative frequencies analysis 
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotation_mtDNA_Clade_Monophyly2.txt -m 5 -p relativeFreq/MtDNA_genes-ML -o results/relativeFreq -g Outgroup
echo "DiscoVista analysis completed. Results are in ${DISCOVISTA_RESULTS}/relativeFreq"



