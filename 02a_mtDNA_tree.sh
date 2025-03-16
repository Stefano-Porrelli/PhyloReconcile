#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#               Script_02a - Generation of Mitochondrial Phylogeny             #
#                                IqTree2                                       #
#------------------------------------------------------------------------------#

# Set paths
BASE_DIR="$(pwd)/PhyloReconcile"
DATA_DIR="${BASE_DIR}/01_initial_data/mtDNA_MSA"
PARTITION="${DATA_DIR}/mtDNA_partitions.txt"
MITO_PHYLO="${BASE_DIR}/02_phylogenies/02a_mtDNA_phylogeny"

# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate PhyloReconcile

# Gather mitochondrial MSA and partitions
mkdir -p "$MITO_PHYLO"
cp "${DATA_DIR}"/mtDNA_MSA.nexus "${MITO_PHYLO}"
cp "${DATA_DIR}"/mtDNA_partitions.txt "${MITO_PHYLO}"
cd "$MITO_PHYLO" || exit

# Infer mtDNA tree with concatenation method
iqtree2 -s mtDNA_MSA.nexus -B 1000 --wbt -nt AUTO --prefix mtDNA_concatenated
# Infer set of mtDNA genes/locus trees
iqtree2 -s mtDNA_MSA.nexus -S mtDNA_partitions.txt -B 1000 --wbt -nt AUTO --prefix mtDNA_loci

# compute gCF and sCF 
iqtree2 -t MtDNA_concatenated.treefile --gcf mtDNA_loci.treefile -s mtDNA_MSA.nexus --scf 100 --prefix MtDNA_concord --cf-verbose --df-tree --cf-quartet
# Cleanup
mkdir -p ./trees
mv ./*.treefile ./trees
mv ./*.tree ./trees
