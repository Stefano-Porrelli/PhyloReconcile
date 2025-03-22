#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#               Script_02a - 1. Generation of Mitochondrial Phylogeny          #
#                                                                              #
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

# Infer mtDNA tree with concatenation method (ML)
iqtree2 -s mtDNA_MSA.nexus -B 1000 --wbt -nt AUTO --prefix MtDNA_concat
cp ./MtDNA_concat.treefile ./MtDNA_concat.tree
# Infer set of mtDNA genes/locus trees
iqtree2 -s mtDNA_MSA.nexus -S mtDNA_partitions.txt -B 1000 --wbt -nt AUTO --prefix MtDNA_loci
# compute gCF and sCF for the concatenated ML mitochondrial tree
iqtree2 -t MtDNA_concat.treefile --gcf MtDNA_loci.treefile -s mtDNA_MSA.nexus --scf 100 --prefix MtDNA_concat_concord --cf-verbose --df-tree --cf-quartet

# ASTRAL: Multispecies coalescent model (MSC)
# Prepare input data files for ASTRAL
ls MtDNA_loci.ufboot > loci_boot.txt

# Run ASTRAL with multilocus bootstrap
astral -i MtDNA_loci.treefile -b loci_boot.txt -r 1000 -o MtDNA_MSC_concord.trees
tail -n 1 MtDNA_MSC_concord.trees > MtDNA_MSC_concord.tree
# compute gCF and sCF 
iqtree2 -t MtDNA_MSC_concord.tree --gcf MtDNA_loci.treefile -s mtDNA_MSA.nexus --scf 100 --prefix MtDNA_MSC_concord --cf-verbose --df-tree --cf-quartet

# Run ASTRAL with local PP
astral -i MtDNA_loci.treefile -o MtDNA_MSC_lpp_concord.tree
# compute gCF and sCF 
iqtree2 -t MtDNA_MSC_lpp_concord.tree --gcf MtDNA_loci.treefile -s mtDNA_MSA.nexus --scf 100 --prefix MtDNA_MSC_lpp_concord --cf-verbose --df-tree --cf-quartet

# Cleanup
mkdir -p ./trees
mv ./*.tree ./trees
rm *.bionj
rm *.ckp.gz
rm *.mldist
rm *.model.gz
rm *.splits.nex
