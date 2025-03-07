#!/bin/bash

BASE_DIR="$(pwd)/Phylogenomic_pipeline"
GENE_TREES_ML="${BASE_DIR}/02_phylogenies/01a_30AX_MLgenetrees"
GENE_TREES_BI="${BASE_DIR}/01_initial_data/Input_files/BI_gene_trees_newick"																				
ASTRAL_DIR="${BASE_DIR}/02_phylogenies/01b_30AX_ML_ASTRAL"
DISCOVISTA_DIR="${BASE_DIR}/05_DiscoVista_analysis"
PARAMETERS="${DISCOVISTA_DIR}/parameters"
DISCOVISTA_TREES="${DISCOVISTA_DIR}/species"
DISCOVISTA_GENES="${DISCOVISTA_DIR}/genetrees"
DISCOVISTA_RESULTS="${DISCOVISTA_DIR}/results"
DISCOVISTA_FREQ="${DISCOVISTA_DIR}/relativeFreq"

# Activate the Conda environment
source activate Phylogenomic_pipeline
# Create necessary directories
mkdir -p "${DISCOVISTA_DIR}" "${PARAMETERS}" "${DISCOVISTA_TREES}" "${DISCOVISTA_RESULTS}"

# Change to the DISCOVISTA_TREES directory
cd "${DISCOVISTA_TREES}" || exit
# Copy all .treefile files from GENE_TREES_ML to DISCOVISTA_TREES for species tree analyses
cp "${GENE_TREES_ML}"/*.treefile "${DISCOVISTA_TREES}"
# Process all .treefile files
for treefile in *.treefile; do
    # Extract the base name of the file (without .treefile extension)
    basename=$(basename "$treefile" .treefile)
    # Extract the gene name 
    gene_name=$(echo "$basename" | cut -d'_' -f2 | cut -d'.' -f1)
    # Create directory name
    dir_name="${gene_name}-ML"
    # Create the directory
    mkdir -p "$dir_name"
    # Move and rename the tree file
    mv "$treefile" "${dir_name}/estimated_species_tree.tree"
    echo "Moved $treefile to ${dir_name}/estimated_species_tree.tree"
done
echo "All trees have been copied, organized into directories, and temporary files cleaned up."

# Copy all .treefile files from GENE_TREES_ML to DISCOVISTA_GENES for gene tree analyses
# Copy and organize ML gene trees
mkdir -p "${DISCOVISTA_GENES}"
cp "${GENE_TREES_ML}"/*.treefile "${DISCOVISTA_GENES}"
# Process all ML .treefile files
cd "${DISCOVISTA_GENES}" || exit
for treefile in "${DISCOVISTA_GENES}"/*.treefile; do
    # Extract the base name of the file (without .treefile extension)
    basename=$(basename "$treefile" .treefile)
    # Extract the gene name 
    gene_name=$(echo "$basename" | cut -d'_' -f2 | cut -d'.' -f1)
    # Create the parent directory for the gene
    mkdir -p "${gene_name}"
    # Create the ML directory within the gene's parent directory
    dir_name="${gene_name}/${gene_name}-ML"
    mkdir -p "$dir_name"
    # Move and rename the tree file
    mv "$treefile" "${dir_name}/estimated_gene_trees.tree"
    echo "Moved $treefile to ${dir_name}/estimated_gene_trees.tree"
done
echo "All trees have been copied, organized into directories, and temporary files cleaned up."

# Copy and organize BI gene trees to DISCOVISTA_GENES for gene tree analyses
cp "${GENE_TREES_BI}"/*.con.tre "${DISCOVISTA_GENES}"											
cd "${DISCOVISTA_GENES}" || exit
for treefile in "${DISCOVISTA_GENES}"/*.con.tre; do
    # Extract the base name of the file (without .con.tre extension)
    basename=$(basename "$treefile" .con.tre)
    # Extract the gene name 
    gene_name=$(echo "$basename" | cut -d'_' -f2 | cut -d'.' -f1)
    # Create the parent directory for the gene
    mkdir -p "${gene_name}"
    # Create the BI directory within the gene's parent directory
    dir_name="${gene_name}/${gene_name}-BI"
    mkdir -p "$dir_name"
    # Move and rename the tree file
    mv "$treefile" "${dir_name}/estimated_gene_trees.tree"
    echo "Moved $treefile to ${dir_name}/estimated_gene_trees.tree"
done
echo "All BI gene trees have been copied and organized into directories."

# Copy and organize all species trees for species tree analyses
cd "${DISCOVISTA_TREES}" || exit
mkdir -p ./ASTRAL-ML
cp "${ASTRAL_DIR}"/30AX_ASTRAL_ML_species.tree ./ASTRAL-ML/estimated_species_tree.tree
# Copy and rename Concatenation and SuperTRI trees											# TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# TODO

# Copy and organize species and gene trees for Relative frequency analysis
cd "${DISCOVISTA_DIR}" || exit
mkdir -p "${DISCOVISTA_FREQ}"/ASTRAL_30genes-ML
# retrieve ASTRAL species tree and ML gene trees
# ML gene trees need to be on a single file. one tree per line
cp "${DISCOVISTA_TREES}"/ASTRAL-ML/estimated_species_tree.tree "${DISCOVISTA_FREQ}"/ASTRAL_30genes-ML
cp "${ASTRAL_DIR}"/ml_best.trees "${DISCOVISTA_FREQ}"/ASTRAL_30genes-ML/estimated_gene_trees.tree

# Pull the DiscoVista Docker image (assume docker is installed)
docker pull esayyari/discovista

# 1 - Discordance analyses on species tree and gene trees - Family/Subfamily level			
# Compare clade support between the different species trees
# and between the individual ML gene trees vs species trees
# 1.1 - Discordance analysis on species tree - Family/Subfamily
# create results directory
cd "${DISCOVISTA_DIR}" || exit
mkdir -p "${DISCOVISTA_RESULTS}"/speciesTree_family_subfamily
# generate clade-defs file for Family/subfamily level (annotation1)
# docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista generate_clade-defs.py parameters/annotation1_FamSubfam.txt parameters/clade-defs1_FamSubfam.txt
# run DiscoVista
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -c parameters/clade-defs1_FamSubfam.txt -p species/ -t 95 -y parameters/newModel.txt -w parameters/newOrders1_FamSubfam.txt -m 0 -o results/speciesTree_family_subfamily
echo "DiscoVista analysis completed. Results are in ${DISCOVISTA_RESULTS}/speciesTree_family_subfamily"
# 1.2 - Discordance analysis on gene trees - Family/Subfamily
# Show the proportion of gene trees for which clades are supported/rejected
mkdir -p "${DISCOVISTA_RESULTS}"/genetrees_family_subfamily
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -c parameters/clade-defs1_FamSubfam.txt -p genetrees/ -t 95 -w parameters/newOrders1_FamSubfam.txt -m 1 -o results/genetrees_family_subfamily
echo "DiscoVista analysis completed. Results are in ${DISCOVISTA_RESULTS}/genetrees_family_subfamily"
# 1.3 - Relative frequencies analysis - Family/Subfamily
# determine frequency of all three topologies around focal branches of the infered species trees
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/Family_subfamily
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotation1_FamSubfam.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/Family_subfamily -g Outgroup
echo "DiscoVista analysis completed. Results are in ${DISCOVISTA_RESULTS}/relativeFreq/Family_subfamily"

# 2 - Discordance analyses on species tree and gene trees - Tribe/Subtribe level
# Compare clade support between the different species trees
# and between the individual ML gene trees vs species trees
# 2.1 - Discordance analysis on species tree - Tribe/Subtribe
# create results directory
cd "${DISCOVISTA_DIR}" || exit
mkdir -p "${DISCOVISTA_RESULTS}"/speciesTree_tribe_subtribe
# generate clade-defs file for Tribe/subtribe level (annotation2)
# docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista generate_clade-defs.py parameters/annotation2_TribeSubtribe.txt parameters/clade-defs2_TribeSubtribe.txt
# run DiscoVista
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -c parameters/clade-defs2_TribeSubtribe.txt -p species/ -t 95 -y parameters/newModel.txt -w parameters/newOrders2_TribeSubtribe.txt -m 0 -o results/speciesTree_tribe_subtribe
echo "DiscoVista analysis completed. Results are in ${DISCOVISTA_RESULTS}/speciesTree_tribe_subtribe"
# 2.2 - Discordance analysis on gene trees - Tribe/Subtribe
# Show the proportion of gene trees for which clades are supported/rejected
mkdir -p "${DISCOVISTA_RESULTS}"/genetrees_tribe_subtribe
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -c parameters/clade-defs2_TribeSubtribe.txt -p genetrees/ -t 95 -w parameters/newOrders2_TribeSubtribe.txt -m 1 -o results/genetrees_tribe_subtribe
echo "DiscoVista analysis completed. Results are in ${DISCOVISTA_RESULTS}/genetrees_tribe_subtribe"
# 1.3 - Relative frequencies analysis - Tribe/subtribe
# determine frequency of all three topologies around focal branches of the infered species trees
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/Tribe_subtribe
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotation2_TribeSubtribe_exclNeotragini.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/Tribe_subtribe -g Outgroup
echo "DiscoVista analysis completed. Results are in ${DISCOVISTA_RESULTS}/relativeFreq/Family_subfamily"

# 3 - Relative frequencies analyses - species level (nodes A-N)
# determine frequency of all three topologies around focal branches 
# that show high level of discordance based on gCF/sCF.
cd "${DISCOVISTA_DIR}" || exit
# 3.1 - Cetacea species level (A)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annA_Cetacea
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationA_Cetacea.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annA_Cetacea -g Outgroup
# 3.2 - Cervidae species level (B)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annB_Cervidae
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationB_Cervidae.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annB_Cervidae -g Outgroup
# 3.3 - Tragelaphini species level (C/D/E)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annC_D_E_Tragelaphini
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationC_D_E_Tragelaphini.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annC_D_E_Tragelaphini -g Outgroup
# 3.4 - Pseudoryna/Bubalina (F)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annF_Pseudoryna
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationF_Pseudoryna.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annF_Pseudoryna -g Outgroup
# 3.5 - Realtionships Bubalina (G)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annG_Bubalina
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationG_Bubalina.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annG_Bubalina -g Outgroup
# 3.6 - Alcelaphini (H)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annH_Alcelaphini
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationH_Alcelaphini.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annH_Alcelaphini -g Outgroup
# 3.7 - Pseudois/Hemitragus/Capra (I)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annI_Pseudois_Hemitragus_Capra
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationI_Pseudois_Hemitragus_Capra.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annI_Pseudois_Hemitragus_Capra -g Outgroup
# 3.8 - intra- and inter-specific relationship within Ovis (J/K)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annJ_K_Ovis
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationJ_K_Ovis.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annJ_K_Ovis -g Outgroup
# 3.9 - Position of Oreotragini (L)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annL_Oreotragini
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationL_Oreotragini.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annL_Oreotragini -g Outgroup
# 3.10 - Position of Oreotragini (M/N)
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/annM_N_Antilopina
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotationM_N_Antilopina.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/annM_N_Antilopina -g Outgroup
# 3.11 - intra- and inter-specific relationship within Bovina
mkdir -p "${DISCOVISTA_RESULTS}"/relativeFreq/ann_Bovina
docker run -v "${DISCOVISTA_DIR}":/data esayyari/discovista discoVista.py -a parameters/annotation_Bovina.txt -m 5 -p relativeFreq/ASTRAL_30genes-ML -o results/relativeFreq/ann_Bovina -g Outgroup


