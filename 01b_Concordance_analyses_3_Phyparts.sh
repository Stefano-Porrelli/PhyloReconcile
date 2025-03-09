#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#                   Script_01b - Concordance Analysis 3                        #
#                Phyparts and quartet-based analyses (ASTRAL)                  #
#------------------------------------------------------------------------------#

PHYPARTS="$HOME/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar"
BASE_DIR="$(pwd)/PhyloReconcile"
CONCORDANCE_DIR="${BASE_DIR}/03_Concordance_analyses"
PHYPARTS_DIR="${CONCORDANCE_DIR}/Phyparts_quartet_analyses"
ASTRAL_DIR="${BASE_DIR}/02_phylogenies/01b_30AX_ML_ASTRAL"
GENE_TREES_ML="${BASE_DIR}/02_phylogenies/01a_30AX_MLgenetrees"
GENE_TREES_ROOTED="${PHYPARTS_DIR}/Gene_trees_rooted"

# Phyparts needs python 2.7 to work, as well as
# matplotlib and ete3
# Create a Conda environment file
cat > "${BASE_DIR}/python27_environment.yml" << EOF
name: python2.7
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=2.7
  - ete3
  - matplotlib
  - newick_utils
EOF

# Create a Conda environment and install required packages
conda env create -y -f "${BASE_DIR}/python27_environment.yml"
# Ensure Conda is initialized properly and activate
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate python2.7

# Install Phyparts
# Needs mvn (Apache Maven) installed on system
# install mvn for MacOS with brew (uncomment next line)
# brew install maven
git clone https://bitbucket.org/blackrim/phyparts.git
cd ./phyparts
./mvn_cmdline.sh

# Create directories and copy files
mkdir -p "${PHYPARTS_DIR}"
mkdir -p "${GENE_TREES_ROOTED}"

cp "${ASTRAL_DIR}/30AX_ASTRAL_ML_species.tree" "${PHYPARTS_DIR}"
cp "${GENE_TREES_ML}"/*.treefile "${GENE_TREES_ROOTED}"

# Loop through all gene trees (.treefile) and root them 
for treefile in "${GENE_TREES_ROOTED}"/*.treefile; do
    # Define output filename
    rooted_tree="${treefile%.treefile}.rooted.tree"
    # Reroot using Newick Utilities (nw_reroot)
    nw_reroot "$treefile" Panthera_leo_Ple1 > "$rooted_tree"
    # Check if the rerooting was successful
    if [ $? -eq 0 ]; then
        echo "Successfully rooted: $treefile -> $rooted_tree"
    else
        echo "Error rerooting: $treefile"
    fi
done
# Remove unrooted trees
rm "${GENE_TREES_ROOTED}"/*.treefile 
cd "${PHYPARTS_DIR}"
# Root species tree
nw_reroot 30AX_ASTRAL_ML_species.tree Panthera_leo_Ple1 > 30AX_ASTRAL_ML_species.rooted.tree

# Run Phyparts
java -jar "$PHYPARTS" -a 1 -v -d "${GENE_TREES_ROOTED}" -m 30AX_ASTRAL_ML_species.rooted.tree -o Phyparts30
# Get pie charts script
curl -O https://raw.githubusercontent.com/mossmatters/MJPythonNotebooks/master/phypartspiecharts.py
# Generate pie charts
python phypartspiecharts.py 30AX_ASTRAL_ML_species.rooted.tree Phyparts30 30

echo "Phyparts completed, species tree with pie charts in pies.svg"