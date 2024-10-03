#!/bin/bash

# Set paths
BASE_DIR="$(pwd)/Phylogenomic_pipeline"
DATA_DIR="${BASE_DIR}/01_initial_data/30AX"
GENE_TREES_DIR="${BASE_DIR}/02_phylogenies/01a_30AX_MLgenetrees"
ASTRAL_DIR="${BASE_DIR}/02_phylogenies/01b_30AX_ASTRAL"
CONCAT_ALIGNMENT="${DATA_DIR}/30AX_concatenated/30AX_concatenated.fasta"
PARTITION="${DATA_DIR}/30AX_concatenated/30AX_partitions.txt"
CONCAT_TREE="${BASE_DIR}/02_phylogenies/01c_30AX_MLconcat"

# Activate the Conda environment
source activate Phylogenomic_pipeline

# 1 - Generate ML gene trees in IQ-TREE 2
mkdir -p "${GENE_TREES_DIR}"
cd "${DATA_DIR}" || exit

for input_file in *.nexus; do
    if [ -f "${input_file}" ]; then
        output_prefix="${GENE_TREES_DIR}/$(basename "${input_file}" .nexus)"
        iqtree2 -s "${input_file}" -B 1000 --wbt -T AUTO -pre "${output_prefix}" || { echo "Error running IQ-TREE 2 for ${input_file}"; exit 1; }
    else
        echo "File not found: ${input_file}"
    fi
done

# 2 - Prepare input data files for ASTRAL
mkdir -p "${ASTRAL_DIR}"
cp "${CONCAT_ALIGNMENT}" "${ASTRAL_DIR}/" || { echo "Error copying concatenated alignment"; exit 1; }

cd "${GENE_TREES_DIR}" || exit
cat *.treefile > "${ASTRAL_DIR}/ml_best.trees" || { echo "Error creating ml_best.trees"; exit 1; }
cp *.ufboot "${ASTRAL_DIR}/" || { echo "Error copying bootstrap files"; exit 1; }

cd "${ASTRAL_DIR}" || exit
ls *.ufboot > ml_boot.txt || { echo "Error creating ml_boot.txt"; exit 1; }

# Check if necessary files exist
if [ ! -f "ml_best.trees" ] || [ ! -f "ml_boot.txt" ] || [ ! -f "${CONCAT_ALIGNMENT}" ]; then
    echo "Error: One or more required files are missing."
    exit 1
fi

# 3 - Run ASTRAL to infer species tree topology and branch lengths (in coalescent units) from the ML gene trees
# bootstrap trees are used to quantify node support with multilocus bootstrapping
astral -i ml_best.trees -b ml_boot.txt -r 1000 -o species_boot.trees || { echo "Error running ASTRAL"; exit 1; }
# The output file "species_boot.trees" contains 1002 lines, 
# the first 1000 lines are the species trees estimated for each of the first 1000 bootstrap trees for each gene 
# line 1001 is a consensus tree
# the last line is the species tree generated from ML gene trees
tail -n 1 species_boot.trees > species.tree || { echo "Error creating species.tree"; exit 1; }

# 4 - Generate gCF and sCF
iqtree2 -t species.tree --gcf ml_best.trees -s "${CONCAT_ALIGNMENT}" --scf 100 --prefix 30AX_ASTRAL_species_tree --cf-verbose --df-tree --cf-quartet || { echo "Error running IQ-TREE 2 for gCF and sCF"; exit 1; }

# Determine quartet support (% of quartets in the gene trees that agree with the branch)
# for the main topology, first alternative and second alternative
astral -i ml_best.trees -b ml_boot.txt -r 1000 -t 8 -o species_boot_t8.trees || { echo "Error running ASTRAL with -t 8"; exit 1; }
# Export branch annotations in .csv (needed for discordance analyses) 
astral -i ml_best.trees -b ml_boot.txt -r 1000 -t 16 -o species_boot_t16.trees || { echo "Error running ASTRAL with -t 8"; exit 1; }

# 5 - Generate ML concatenated phylogeny IQ-TREE 2
mkdir -p "${CONCAT_TREE}"
cd "${CONCAT_TREE}" || exit
iqtree2 -s "${CONCAT_ALIGNMENT}" -p "${PARTITION}" --prefix 30AX_MLconcatenation -B 1000 -T AUTO || { echo "Error running IQ-TREE 2 for ML concatenation"; exit 1; }

echo "ASTRAL species tree with bootstrap and Concordance values for node support generated."
echo "ML concatenation phylogeny generated."
