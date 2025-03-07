#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#           Script_01a - Generation of gene trees and species tree under       #
#                the multi-species coalescent (MSC) model in ASTRAL            #
#------------------------------------------------------------------------------#

# Set paths
BASE_DIR="$(pwd)/PhyloReconcile"
DATA_DIR="${BASE_DIR}/01_initial_data/30AX_MSAs"
ML_GENE_TREES_DIR="${BASE_DIR}/02_phylogenies/01a_30AX_MLgenetrees"
ML_ASTRAL_DIR="${BASE_DIR}/02_phylogenies/01b_30AX_ML_ASTRAL"
BI_GENE_TREES_DIR="${BASE_DIR}/02_phylogenies/01c_30AX_BIgenetrees"
BI_ASTRAL_DIR="${BASE_DIR}/02_phylogenies/01d_30AX_BI_ASTRAL"
CONCAT_ALIGNMENT="${DATA_DIR}/30AX_concatenated/30AX_concatenated.fasta"
PARTITION="${DATA_DIR}/30AX_concatenated/30AX_partitions.txt"
CONCAT_TREE="${BASE_DIR}/02_phylogenies/01e_30AX_MLconcat"
MRBAYES_BLOCKS="${DATA_DIR}/30AX_MrBayes_blocks"

# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate PhyloReconcile

# 1 - ASTRAL MSC ML Phylogeny
# 1.1 - Generate ML gene trees in IQ-TREE 2
cd "${DATA_DIR}" || exit

for input_file in *.nexus; do
    if [ -f "${input_file}" ]; then
        output_prefix="${ML_GENE_TREES_DIR}/$(basename "${input_file}" .nexus)"
        iqtree2 -s "${input_file}" -B 1000 --wbt -T AUTO -pre "${output_prefix}" || { echo "Error running IQ-TREE 2 for ${input_file}"; exit 1; }
    else
        echo "File not found: ${input_file}"
    fi
done

# 1.2 - Prepare input data files for ASTRAL
cp "${CONCAT_ALIGNMENT}" "${ML_ASTRAL_DIR}/" || { echo "Error copying concatenated alignment"; exit 1; }
cd "${ML_GENE_TREES_DIR}" || exit
cat *.treefile > "${ML_ASTRAL_DIR}/ml_best.trees" || { echo "Error creating ml_best.trees"; exit 1; }
cp *.ufboot "${ML_ASTRAL_DIR}/" || { echo "Error copying bootstrap files"; exit 1; }
cd "${ML_ASTRAL_DIR}" || exit
ls *.ufboot > ml_boot.txt || { echo "Error creating ml_boot.txt"; exit 1; }

# Check if necessary files exist
if [ ! -f "ml_best.trees" ] || [ ! -f "ml_boot.txt" ] || [ ! -f "${CONCAT_ALIGNMENT}" ]; then
    echo "Error: One or more required files are missing."
    exit 1
fi

# 1.3 - Run ASTRAL to infer species tree topology and branch lengths (in coalescent units) from the ML gene trees
astral -i ml_best.trees -b ml_boot.txt -r 1000 -o species_boot.trees || { echo "Error running ASTRAL"; exit 1; }
# the last line of the output file "species_boot.trees" is the species tree generated from ML gene trees
tail -n 1 species_boot.trees > 30AX_ASTRAL_ML_species.tree || { echo "Error creating species.tree"; exit 1; }

# 1.4 - Generate gCF and sCF
iqtree2 -t 30AX_ASTRAL_ML_species.tree --gcf ml_best.trees -s "${CONCAT_ALIGNMENT}" --scf 100 --prefix 30AX_ASTRAL_ML_species_tree --cf-verbose --df-tree --cf-quartet || { echo "Error running IQ-TREE 2 for gCF and sCF"; exit 1; }
# Determine quartet support (% of quartets in the gene trees that agree with the branch)
# for the main topology, first alternative and second alternative
astral -i ml_best.trees -b ml_boot.txt -r 1000 -t 8 -o species_boot_t8.trees || { echo "Error running ASTRAL with -t 8"; exit 1; }
# Export branch annotations in tab delimited file
astral -i ml_best.trees -b ml_boot.txt -r 1000 -t 16 -o species_boot_t16.trees || { echo "Error running ASTRAL with -t 16"; exit 1; }

echo "ASTRAL species tree with bootstrap and Concordance values for node support generated."

# 2 - ASTRAL MSC BI Phylogeny
# 2.1 - Generate BI gene trees in MrBayes 
cd "${BI_GENE_TREES_DIR}" || exit
cp "${MRBAYES_BLOCKS}"/*.nexus "${BI_GENE_TREES_DIR}"

for input_file in *.nexus; do
    if [ -f "${input_file}" ]; then
        output_prefix="$(basename "${input_file}" .nexus)"
        mb "${input_file}" || { echo "Error running MrBayes for ${input_file}"; exit 1; }
    else
        echo "File not found: ${input_file}"
    fi
done

# 2.2 - Prepare input data files for ASTRAL
cp "${CONCAT_ALIGNMENT}" "${BI_ASTRAL_DIR}/" || { echo "Error copying concatenated alignment"; exit 1; }
# convert nexus to newick with AfterPhylo
curl -O https://raw.githubusercontent.com/qiyunzhu/AfterPhylo/master/AfterPhylo.pl
chmod +x AfterPhylo.pl

for treefile in *.con.tre; do
    perl AfterPhylo.pl -confonly -format=newick "$treefile"
done

cp "${BI_GENE_TREES_DIR}"/*out.tre "${BI_ASTRAL_DIR}" 
cd "${BI_ASTRAL_DIR}"

for treefile in *.tre; do
    # Ensure there's a newline at the end of the file, then append it to bi.trees
    echo "" >> "$treefile"
    cat "$treefile" >> "${BI_ASTRAL_DIR}/bi.trees"
done || { echo "Error creating bi.trees"; exit 1; }

# 2.3 - Run ASTRAL to infer species tree topology and branch lengths (in coalescent units) from the ML gene trees
astral -i bi.trees -r 1000 -o 30AX_ASTRAL_BI_species.tree  || { echo "Error running ASTRAL"; exit 1; }
# 2.4 - Generate gCF and sCF
iqtree2 -t 30AX_ASTRAL_BI_species.tree --gcf bi.trees -s "${CONCAT_ALIGNMENT}" --scf 100 --prefix 30AX_ASTRAL_BI_species_tree_gCF --cf-verbose --df-tree --cf-quartet || { echo "Error running IQ-TREE 2 for gCF and sCF"; exit 1; }
# Determine quartet support (% of quartets in the gene trees that agree with the branch)
# for the main topology, first alternative and second alternative
astral -i bi.trees -r 1000 -t 8 -o species_bi_t8.trees || { echo "Error running ASTRAL with -t 8"; exit 1; }
# Export branch annotations in tab delimited file
astral -i bi.trees -r 1000 -t 16 -o species_bi_t16.trees || { echo "Error running ASTRAL with -t 16"; exit 1; }

echo "ASTRAL species tree with PP and Concordance values for node support generated."

# 3 - Generate ML concatenated phylogeny IQ-TREE 2
cd "${CONCAT_TREE}" || exit
iqtree2 -s "${CONCAT_ALIGNMENT}" -p "${PARTITION}" --prefix 30AX_ML_concatenation -B 1000 -T AUTO || { echo "Error running IQ-TREE 2 for ML concatenation"; exit 1; }

echo "ML concatenation phylogeny generated."
