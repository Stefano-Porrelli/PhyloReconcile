#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#           Script_00 - Installation of environment and packages,              #
#                   directories structure and data retrieval                   #
#------------------------------------------------------------------------------#

# Set base directory
BASE_DIR="$(pwd)/PhyloReconcile"
# Make directories and subdirectories
mkdir -p "${BASE_DIR}"
mkdir -p \
  "${BASE_DIR}/02_phylogenies/01a_30AX_MLgenetrees" \
  "${BASE_DIR}/02_phylogenies/01b_30AX_ML_ASTRAL" \
  "${BASE_DIR}/02_phylogenies/01c_30AX_BIgenetrees" \
  "${BASE_DIR}/02_phylogenies/01d_30AX_BI_ASTRAL" \
  "${BASE_DIR}/02_phylogenies/01e_30AX_MLconcat" \
  "${BASE_DIR}/03_Concordance_analyses/Concordance_factors"

# Install Miniconda: uncomment next block for installation
# mkdir -p ~/miniconda3
# curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
# bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
# rm ~/miniconda3/miniconda.sh
# source ~/miniconda3/bin/activate
# conda init --all

# Create a Conda environment file
cat > "${BASE_DIR}/environment.yml" << EOF
name: PhyloReconcile
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
conda activate PhyloReconcile

# Install snp-sites, samtools, vcftools and bcftools,
# newick_utils, and numpy
conda install -y bioconda::snp-sites
conda install -y bioconda::samtools
conda install -y bioconda::vcftools
conda install -y soil::bcftools
conda install -y bioconda::newick_utils
conda install -y conda-forge::numpy

# Retrieve data
git clone --depth 1 https://github.com/Stefano-Porrelli/PhyloReconcile.git temp_repo
cp -r temp_repo/Datasets/* "${DATA_DIR}/"
rm -rf temp_repo

echo "Setup and data preparation completed."
