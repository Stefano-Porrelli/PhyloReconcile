#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#                   Script_02b - 1. Detection of ILS with PhyParts             #
#                       Mitochondrial and Y_chromosome data                    #
#------------------------------------------------------------------------------#

PHYPARTS="$HOME/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar"
BASE_DIR="$(pwd)/PhyloReconcile"
CONCORDANCE_DIR="${BASE_DIR}/07_Phyparts_MtDNA_MSY"

MTDNA_PHYPARTS="${CONCORDANCE_DIR}/MtDNA"
MTDNA_PHYLO="${BASE_DIR}/02_phylogenies/02a_mtDNA_phylogeny/trees/MtDNA_concat.tree"
GENE_TREES_MTDNA="${BASE_DIR}/02_phylogenies/02a_mtDNA_phylogeny/MtDNA_loci.treefile"
GENE_TREES_ROOTED_MTDNA="${MTDNA_PHYPARTS}/Gene_trees_rooted"

MSY_PHYPARTS="${CONCORDANCE_DIR}/MSYDNA"
MSY_PHYLO="${BASE_DIR}/02_phylogenies/02b_MSY_phylogeny/trees/MSY_concatenation.tree"
GENE_TREES_MSY="${BASE_DIR}/02_phylogenies/02b_MSY_phylogeny/Yloci.treefile"
GENE_TREES_ROOTED_MSY="${MSY_PHYPARTS}/Gene_trees_rooted"

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

# Mitochondrial data
# Create directories and copy files
mkdir -p "${MTDNA_PHYPARTS}"
mkdir -p "${GENE_TREES_ROOTED_MTDNA}"
cp "${MTDNA_PHYLO}" "${MTDNA_PHYPARTS}"
cp "${GENE_TREES_MTDNA}" "${GENE_TREES_ROOTED_MTDNA}"
# Loop through all gene trees (.treefile) and root them 
for treefile in "${GENE_TREES_ROOTED_MTDNA}"/*.treefile; do
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
rm "${GENE_TREES_ROOTED_MTDNA}"/*.treefile 
cd "${MTDNA_PHYPARTS}"
# Root species tree
nw_reroot MtDNA_concat.tree Panthera_leo_Ple1 > MtDNA_concat.rooted.tree
# Run Phyparts
java -jar "$PHYPARTS" -a 1 -v -d "${GENE_TREES_ROOTED_MTDNA}" -m MtDNA_concat.rooted.tree -o Phyparts14
# Get pie charts script
curl -O https://raw.githubusercontent.com/mossmatters/MJPythonNotebooks/master/phypartspiecharts.py
# Generate pie charts
python phypartspiecharts.py MtDNA_concat.rooted.tree Phyparts14 14
echo "Phyparts completed, species tree with pie charts in pies.svg"


# Y chromosome data
# Create directories and copy files
mkdir -p "${MSY_PHYPARTS}"
mkdir -p "${GENE_TREES_ROOTED_MSY}"
cp "${MSY_PHYLO}" "${MSY_PHYPARTS}"
cp "${GENE_TREES_MSY}" "${GENE_TREES_ROOTED_MSY}"
# Loop through all gene trees (.treefile) and root them 
for treefile in "${GENE_TREES_ROOTED_MSY}"/*.treefile; do
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
rm "${GENE_TREES_ROOTED_MSY}"/*.treefile 
cd "${MSY_PHYPARTS}"
# Root species tree
nw_reroot MSY_concatenation.tree Panthera_leo_Ple1 > MSY_concatenation.rooted.tree
# Run Phyparts
java -jar "$PHYPARTS" -a 1 -v -d "${GENE_TREES_ROOTED_MSY}" -m MSY_concatenation.rooted.tree -o Phyparts2
# Get pie charts script
curl -O https://raw.githubusercontent.com/mossmatters/MJPythonNotebooks/master/phypartspiecharts.py
# Generate pie charts
python phypartspiecharts.py MSY_concatenation.rooted.tree Phyparts2 2
echo "Phyparts completed, species tree with pie charts in pies.svg"
