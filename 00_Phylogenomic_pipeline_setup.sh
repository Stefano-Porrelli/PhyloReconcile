#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#           Script_00 - Installation of environment and packages,              #
#                   directories structure and data retrieval                   #
#------------------------------------------------------------------------------#

# Set base directory
BASE_DIR="$(pwd)/Phylogenomic_pipeline"
# Make directories and subdirectories
mkdir -p "${BASE_DIR}"
mkdir -p \
  "${BASE_DIR}/02_phylogenies/01a_30AX_MLgenetrees" \
  "${BASE_DIR}/02_phylogenies/01b_30AX_ML_ASTRAL" \
  "${BASE_DIR}/02_phylogenies/01c_30AX_BIgenetrees" \
  "${BASE_DIR}/02_phylogenies/01d_30AX_BI_ASTRAL" \
  "${BASE_DIR}/02_phylogenies/01e_30AX_MLconcat" \
  "${BASE_DIR}/R_scripts" \
  "${BASE_DIR}/03_Concordance_analyses/Concordance_factors"

# Install Miniconda
# (Instructions for installing Miniconda go here)

# Create a Conda environment file
cat > "${BASE_DIR}/environment.yml" << EOF
name: Phylogenomic_pipeline
channels:
  - conda-forge
  - bioconda
dependencies:
  - python
  - iqtree
  - astral-tree
  - mrbayes
  - amas
  - r-base
  - r-viridis
  - r-ggplot2
  - r-dplyr
  - r-ggrepel
  - r-ggally
  - r-entropy
  - r-patchwork
EOF

# Create a Conda environment and install required packages
conda env create -f "${BASE_DIR}/environment.yml"
# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate Phylogenomic_pipeline


# Install snp-sites, samtools, vcftools and bcftools,
# needed for Dsuite analysis, comment next block if not needed
conda install -y bioconda::snp-sites
conda install -y bioconda::samtools
conda install -y bioconda::vcftools
conda install -y soil::bcftools
conda install -y bioconda::newick_utils

# Retrieve data
git clone --depth 1 https://github.com/Stefano-Porrelli/Phylogenomic_pipeline.git temp_repo
cp -r temp_repo/Datasets/* "${DATA_DIR}/"
rm -rf temp_repo

echo "Setup and data preparation completed."
