#!/bin/bash

# Set base directory
BASE_DIR="$(pwd)/Phylogenomic_pipeline"
DATA_DIR="${BASE_DIR}/01_initial_data"
CONCAT_DIR="${DATA_DIR}/30AX_concatenated"

# Make directories and subdirectories
mkdir -p "${BASE_DIR}" "${DATA_DIR}"
mkdir -p "${BASE_DIR}/02_phylogenies/01a_30AX_MLgenetrees" "${BASE_DIR}/02_phylogenies/01b_30AX_ASTRAL" "${BASE_DIR}/02_phylogenies/01c_30AX_MLconcat" "${BASE_DIR}/R_scripts" "${BASE_DIR}/03_CF_DF_analysis/CF_plots"

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
  - amas
  - r-base
  - r-viridis
  - r-ggplot2
  - r-dplyr
  - r-ggrepel
  - r-ggally
  - r-entropy
EOF

# Create a Conda environment and install required packages
conda env create -f "${BASE_DIR}/environment.yml"
source activate Phylogenomic_pipeline

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
