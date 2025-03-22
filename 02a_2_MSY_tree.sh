#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#          Script_022 - 2. Generation of MSY Phylogeny (AMELY/DDX3Y)           #
#                                                                              #
#------------------------------------------------------------------------------#

# Set paths
BASE_DIR="$(pwd)/PhyloReconcile"
DATA_DIR="${BASE_DIR}/01_initial_data/YDNA_MSA"
PARTITION="${DATA_DIR}/MSY_partitions.txt"
MSY_PHYLO="${BASE_DIR}/02_phylogenies/02b_MSY_phylogeny"

# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate PhyloReconcile

# Gather data: Y chromosome and X-homologs MSA and partitions
mkdir -p "$MSY_PHYLO"
cp "${DATA_DIR}"/*.nexus "${MSY_PHYLO}"
cp "${PARTITION}" "${MSY_PHYLO}"
cd "$MSY_PHYLO" || exit

# Compare X and Y markers in a phylogenetic context
# to validate Y chromosome sequence authenticity
mkdir ./XY_comparison
mv 31_AMEL_MSA.nexus ./XY_comparison
mv 32_DDX3_MSA.nexus ./XY_comparison
cd XY_comparison
# Infer tree for the combined X and Y sequences of AMEL
iqtree2 -s 31_AMEL_MSA.nexus -B 1000 --wbt -nt AUTO --prefix AMEL_XY_comparison
# Infer tree for the combined X and Y sequences of DDX3
iqtree2 -s 32_DDX3_MSA.nexus -B 1000 --wbt -nt AUTO --prefix DDX3_XY_comparison
# Cleanup
mkdir -p ./trees
mv ./AMEL_XY_comparison.treefile ./trees/AMEL_XY_comparison.tree
mv ./DDX3_XY_comparison.treefile ./trees/DDX3_XY_comparison.tree
rm *.bionj
rm *.ckp.gz
rm *.mldist
rm *.model.gz
rm *.splits.nex

# MSY Dataset
# Generate concatenated reference tree
cd "${MSY_PHYLO}"
iqtree2 -s 33_MSY_MSA.nexus -p MSY_partitions.txt --prefix MSY_concatenation -B 1000 --wbt -nt AUTO
# Generate individual loci trees
iqtree2 -s 33_MSY_MSA.nexus -S MSY_partitions.txt -B 1000 --wbt --prefix Yloci -nt AUTO
# compute gCF and sCF for the concatenated ML mitochondrial tree
iqtree2 -t MSY_concatenation.treefile --gcf Yloci.treefile -s 33_MSY_MSA.nexus --scf 100 --prefix MSY_concat_concord --cf-verbose --df-tree --cf-quartet

# ASTRAL: Multispecies coalescent model (MSC)
# Prepare input data files for ASTRAL
ls Yloci.ufboot > Yloci_boot.txt

# Run ASTRAL with multilocus bootstrap
astral -i Yloci.treefile -b Yloci_boot.txt -r 1000 -o MSY_MSC_concord.trees
tail -n 1 MSY_MSC_concord.trees > MSY_MSC_concord.tree
# compute gCF and sCF 
iqtree2 -t MSY_MSC_concord.tree --gcf Yloci.treefile -s 33_MSY_MSA.nexus --scf 100 --prefix MSY_MSC_concord --cf-verbose --df-tree --cf-quartet

# Run ASTRAL with local PP
astral -i Yloci.treefile -o MSY_MSC_lpp_concord.tree
# compute gCF and sCF 
iqtree2 -t MSY_MSC_lpp_concord.tree --gcf Yloci.treefile -s 33_MSY_MSA.nexus --scf 100 --prefix MSY_MSC_lpp_concord --cf-verbose --df-tree --cf-quartet

# Cleanup
mkdir -p ./trees
mv ./*.tree ./trees
mv ./MSY_concatenation.treefile ./trees/MSY_concatenation.tree
rm *.bionj
rm *.ckp.gz
rm *.mldist
rm *.model.gz
rm *.splits.nex

