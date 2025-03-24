#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#          Script_02a - 3. comparison of topological discordance               #
#                         between inheritance modes                            #
#------------------------------------------------------------------------------#

# Set paths
BASE_DIR="$(pwd)/PhyloReconcile"
YLOCI="${BASE_DIR}/02_phylogenies/02b_MSY_phylogeny/Yloci.treefile"
MSY="${BASE_DIR}/02_phylogenies/02b_MSY_phylogeny/trees/MSY_concatenation.tree"
MTDNA="${BASE_DIR}/02_phylogenies/02a_mtDNA_phylogeny/trees/MtDNA_concat.tree"
NUDNA="${BASE_DIR}/02_phylogenies/01b_30AX_ML_ASTRAL/30AX_ASTRAL_ML_species.tree"
RES_DIR="${BASE_DIR}/06_comparison_inheritance_modes"

# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate PhyloReconcile

# Gather data: Y chromosome and X-homologs MSA and partitions
mkdir -p "$RES_DIR"
cd "$RES_DIR" || exit
cp "${YLOCI}" "${RES_DIR}"
cp "${MTDNA}" "${RES_DIR}"
cp "${NUDNA}" "${RES_DIR}"
cp "${MSY}" "${RES_DIR}"
sed -n '1p' Yloci.treefile > AMELY.tree
sed -n '2p' Yloci.treefile > DDX3Y.tree
rm Yloci.treefile

# 1. Compare AMELY and DDX3Y with Cophylo
# Create a separate R script
cat > 01_AMELY_DDX3Y_tanglegram.R << 'EOF'
# Load required libraries
library(phytools)
library(ape)
# Read the trees
AMELY <- read.tree("AMELY.tree")
DDX3Y <- read.tree("DDX3Y.tree")
# Use reorder.phylo to sort by node values
AMELY <- reorder.phylo(AMELY, "postorder")
DDX3Y <- reorder.phylo(DDX3Y, "postorder")
# Apply ladderize with right=TRUE for descending order
AMELY <- ladderize(AMELY, right = TRUE)
DDX3Y <- ladderize(DDX3Y, right = TRUE)
# Remove branch lengths (uncomment)
# AMELY <- compute.brlen(AMELY)
# DDX3Y <- compute.brlen(DDX3Y)
# Create Cophylo object
obj <- cophylo(AMELY, DDX3Y, rotate = TRUE)

# Open PDF device with dimensions 14.2 x 8.27 inches
pdf("AMELY_DDX3Y_comparison.pdf", width=14.2, height=16.27)

# Plot the trees
plot(obj, type="phylogram", fsize=0.5, ftype="i", part=0.48, lwd=0.5,
     mar=c(0.1,0.1,2.1,0.1),
     tip.len=0.001, tip.lwd=0.35,
     link.type="curved", link.lty="solid", link.lwd=1,
     link.col=make.transparent("green",0.25), pts=FALSE)
# Add labels to the plot
mtext("a) AMELY", at=-0.5, adj=0)
mtext("b) DDX3Y", at=0.43, adj=0)
# Report node support below 100 for left tree (AMELY)
left_node_labels <- as.numeric(obj$trees[[1]]$node.label[2:Nnode(obj$trees[[1]])])
left_node_numbers <- 2:Nnode(obj$trees[[1]]) + Ntip(obj$trees[[1]])
left_show_labels <- left_node_labels < 100
left_rounded_labels <- round(left_node_labels[left_show_labels], 0)
nodelabels.cophylo(left_rounded_labels,
                   left_node_numbers[left_show_labels],
                   frame="none",
                   cex=0.4,
                   adj=c(1,-0.4),
                   which="left")
# Report node support below 100 for right tree (DDX3Y)
right_node_labels <- as.numeric(obj$trees[[2]]$node.label[2:Nnode(obj$trees[[2]])])
right_node_numbers <- 2:Nnode(obj$trees[[2]]) + Ntip(obj$trees[[2]])
right_show_labels <- right_node_labels < 100
right_rounded_labels <- round(right_node_labels[right_show_labels], 0)
nodelabels.cophylo(right_rounded_labels,
                   right_node_numbers[right_show_labels],
                   frame="none",
                   cex=0.4,
                   adj=c(0,-0.4),
                   which="right")
dev.off()
EOF
# Execute the R script
Rscript 01_AMELY_DDX3Y_tanglegram.R
echo "Analysis complete. Output saved to ${RES_DIR}/AMELY_DDX3Y_comparison.pdf"

# 2. Compare mitochondrial and nuclear phylogenies with Cophylo
# Create a separate R script
cat > 02_MitoNuclear_tanglegram.R << 'EOF'
# Load required libraries
library(phytools)
library(ape)
# read trees
nuDNA <- read.tree("30AX_ASTRAL_ML_species.tree")
mtDNA <- read.tree("MtDNA_concat.tree")
# Use reorder.phylo to sort by node values
nuDNA <- reorder.phylo(nuDNA, "postorder")
mtDNA <- reorder.phylo(mtDNA, "postorder")
# Apply ladderize with right=TRUE for descending order
nuDNA <- ladderize(nuDNA, right = TRUE)
mtDNA <- ladderize(mtDNA, right = TRUE)
# Remove branch lengths (uncomment)
nuDNA <- compute.brlen(nuDNA)
mtDNA <- compute.brlen(mtDNA)
# Create Cophylo object
obj <- cophylo(mtDNA,nuDNA, rotate = TRUE)

# Open PDF device with dimensions 14.2 x 8.27 inches
pdf("mtDNA_nuDNA_comparison.pdf", width=14.2, height=16.27)

# Plot the trees
plot(obj, type="phylogram", fsize=0.5, ftype="i", part=0.48, lwd=0.5,
     mar=c(0.1,0.1,2.1,0.1),
     tip.len=0.001, tip.lwd=0.35,
     link.type="curved", link.lty="solid", link.lwd=1,
     link.col=make.transparent("green",0.5), pts=FALSE)
# Report node support below 100 for left tree
left_node_labels <- as.numeric(obj$trees[[1]]$node.label[2:Nnode(obj$trees[[1]])])
left_node_numbers <- 2:Nnode(obj$trees[[1]]) + Ntip(obj$trees[[1]])
# Ensure we're strictly filtering values that are equal to or greater than 100
left_show_labels <- left_node_labels < 100
left_rounded_labels <- left_node_labels[left_show_labels]
nodelabels.cophylo(left_rounded_labels,
                  left_node_numbers[left_show_labels],
                  frame="none",
                  cex=0.4,
                  adj=c(1,-0.4),
                  which="left")
# Report node support below 100 for right tree
right_node_labels <- as.numeric(obj$trees[[2]]$node.label[2:Nnode(obj$trees[[2]])])
right_node_numbers <- 2:Nnode(obj$trees[[2]]) + Ntip(obj$trees[[2]])
# Ensure we're strictly filtering values that are equal to or greater than 100
right_show_labels <- right_node_labels < 100
right_rounded_labels <- right_node_labels[right_show_labels]
nodelabels.cophylo(right_rounded_labels,
                  right_node_numbers[right_show_labels],
                  frame="none",
                  cex=0.4,
                  adj=c(0,-0.4),
                  which="right")

# Add labels to the plot
mtext("a) mtDNA", at=-0.5, adj=0)
mtext("b) nuDNA", at=0.43, adj=0)
dev.off()
EOF

# Execute the R script
Rscript 02_MitoNuclear_tanglegram.R
echo "Analysis complete. Output saved to ${RES_DIR}/mtDNA_nuDNA_comparison.pdf"

# 3. Compare mitochondrial and MSY phylogenies with Cophylo
# Create a separate R script
cat > 03_MitoNuclearMSY_tanglegram.R << 'EOF'
# Load required libraries
library(phytools)
library(ape)
# read trees
mtDNA <- read.tree("MtDNA_concat.tree")
msyDNA <- read.tree("MSY_concatenation.tree")
# Use reorder.phylo to sort by node values
mtDNA <- reorder.phylo(mtDNA, "postorder")
msyDNA <- reorder.phylo(msyDNA, "postorder")
# Apply ladderize with right=TRUE for descending order
mtDNA <- ladderize(mtDNA, right = TRUE)
msyDNA <- ladderize(msyDNA, right = TRUE)
# Remove branch lengths (uncomment)
mtDNA <- compute.brlen(mtDNA)
msyDNA <- compute.brlen(msyDNA)
# Create Cophylo object
obj <- cophylo(mtDNA,msyDNA, rotate = TRUE)

# Open PDF device with dimensions 14.2 x 8.27 inches
pdf("mtDNA_msyDNA_comparison.pdf", width=14.2, height=16.27)

# Plot the trees
plot(obj, type="phylogram", fsize=0.5, ftype="i", part=0.48, lwd=0.5,
     mar=c(0.1,0.1,2.1,0.1),
     tip.len=0.001, tip.lwd=0.35,
     link.type="curved", link.lty="solid", link.lwd=1,
     link.col=make.transparent("green",0.5), pts=FALSE)
# Report node support below 100 for left tree
left_node_labels <- as.numeric(obj$trees[[1]]$node.label[2:Nnode(obj$trees[[1]])])
left_node_numbers <- 2:Nnode(obj$trees[[1]]) + Ntip(obj$trees[[1]])
# Ensure we're strictly filtering values that are equal to or greater than 100
left_show_labels <- left_node_labels < 100
left_rounded_labels <- left_node_labels[left_show_labels]
nodelabels.cophylo(left_rounded_labels,
                  left_node_numbers[left_show_labels],
                  frame="none",
                  cex=0.4,
                  adj=c(1,-0.4),
                  which="left")
# Report node support below 100 for right tree
right_node_labels <- as.numeric(obj$trees[[2]]$node.label[2:Nnode(obj$trees[[2]])])
right_node_numbers <- 2:Nnode(obj$trees[[2]]) + Ntip(obj$trees[[2]])
# Ensure we're strictly filtering values that are equal to or greater than 100
right_show_labels <- right_node_labels < 100
right_rounded_labels <- right_node_labels[right_show_labels]
nodelabels.cophylo(right_rounded_labels,
                  right_node_numbers[right_show_labels],
                  frame="none",
                  cex=0.4,
                  adj=c(0,-0.4),
                  which="right")
# Add labels to the plot
mtext("a) mtDNA", at=-0.5, adj=0)
mtext("b) msyDNA", at=0.43, adj=0)
dev.off()
EOF

# Execute the R script
Rscript 03_MitoNuclearMSY_tanglegram.R
echo "Analysis complete. Output saved to ${RES_DIR}/mtDNA_msyDNA_comparison.pdf"

# 4. Compare MSY and nuDNA phylogenies with Cophylo
# Create a separate R script
cat > 04_MSYDNA_nuDNA_tanglegram.R << 'EOF'
# Load required libraries
library(phytools)
library(ape)
# read trees
msyDNA <- read.tree("MSY_concatenation.tree")
nuDNA <- read.tree("30AX_ASTRAL_ML_species.tree")
# Use reorder.phylo to sort by node values
msyDNA <- reorder.phylo(msyDNA, "postorder")
nuDNA <- reorder.phylo(nuDNA, "postorder")
# Apply ladderize with right=TRUE for descending order
msyDNA <- ladderize(msyDNA, right = TRUE)
nuDNA <- ladderize(nuDNA, right = TRUE)
# Remove branch lengths (uncomment)
msyDNA <- compute.brlen(msyDNA)
nuDNA <- compute.brlen(nuDNA)
# Create Cophylo object
obj <- cophylo(msyDNA,nuDNA, rotate = TRUE)

# Open PDF device with dimensions 14.2 x 8.27 inches
pdf("msyDNA_nuDNA_comparison.pdf", width=14.2, height=16.27)

# Plot the trees
plot(obj, type="phylogram", fsize=0.5, ftype="i", part=0.48, lwd=0.5,
     mar=c(0.1,0.1,2.1,0.1),
     tip.len=0.001, tip.lwd=0.35,
     link.type="curved", link.lty="solid", link.lwd=1,
     link.col=make.transparent("green",0.5), pts=FALSE)
# Report node support below 100 for left tree
left_node_labels <- as.numeric(obj$trees[[1]]$node.label[2:Nnode(obj$trees[[1]])])
left_node_numbers <- 2:Nnode(obj$trees[[1]]) + Ntip(obj$trees[[1]])
# Ensure we're strictly filtering values that are equal to or greater than 100
left_show_labels <- left_node_labels < 100
left_rounded_labels <- left_node_labels[left_show_labels]  # No rounding as requested
nodelabels.cophylo(left_rounded_labels,
                  left_node_numbers[left_show_labels],
                  frame="none",
                  cex=0.4,
                  adj=c(1,-0.4),
                  which="left")
# Report node support below 100 for right tree
right_node_labels <- as.numeric(obj$trees[[2]]$node.label[2:Nnode(obj$trees[[2]])])
right_node_numbers <- 2:Nnode(obj$trees[[2]]) + Ntip(obj$trees[[2]])
# Ensure we're strictly filtering values that are equal to or greater than 100
right_show_labels <- right_node_labels < 100
right_rounded_labels <- right_node_labels[right_show_labels]  # No rounding as requested
nodelabels.cophylo(right_rounded_labels,
                  right_node_numbers[right_show_labels],
                  frame="none",
                  cex=0.4,
                  adj=c(0,-0.4),
                  which="right")
# Add labels to the plot
mtext("a) msyDNA", at=-0.5, adj=0)
mtext("b) nuDNA", at=0.43, adj=0)
dev.off()
EOF

# Execute the R script
Rscript 04_MSYDNA_nuDNA_tanglegram.R
echo "Analysis complete. Output saved to ${RES_DIR}/msyDNA_nuDNA_comparison.pdf"
