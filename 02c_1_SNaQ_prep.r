#!/usr/bin/env Rscript
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#                      Script_02c - 1. Data prep for SNaQ                      #
#                                                                              # 
# the following inputs files are needed (examples in 1_initial data):          #
# [group]_individuals.txt = list of individuals in clade                       #
# [group]_species.txt = list of species in clade                               #
# [group]_map.csv = mapping of individuals to species                          #
#                                                                              #
# All other inputs (gene trees, species tree) are generated in previous scripts#
#------------------------------------------------------------------------------#

# 1. Load required libraries
library('MSCquartets')
library('ape')
library('gtools')
library('stringr')

# 2. Set the group to analyze
group = " " # Name of the group (based on input files and directory, eg. Tragelaphini)

# 3. Setup directory structure
BASE_DIR <- file.path(getwd(), "PhyloReconcile")
NUDNA_PHYLO <- file.path(BASE_DIR, "02_phylogenies", "01b_30AX_ML_ASTRAL", "30AX_ASTRAL_ML_species.tree")
GENETREES_NUDNA <- file.path(BASE_DIR, "02_phylogenies", "01b_30AX_ML_ASTRAL", "ml_best.trees")
INPUTS <- file.path(BASE_DIR, "01_initial_data", "input_files", "SNaQ_NANUQ")
MSC_RESULTS <- file.path(BASE_DIR, "09_Network_Analyses")

# 4. Create results directory based on group name
group_dir <- file.path(MSC_RESULTS, group)
if (!dir.exists(group_dir)) {
  dir.create(group_dir, recursive = TRUE)
}

# 5. Load input files
# Move to the input directory where our text files are located
setwd(INPUTS)
species <- read.delim(paste0(group, "_species.txt"), header = FALSE)
indivs <- read.delim(paste0(group, "_individuals.txt"), header = FALSE)
outgroup <- species[length(species)]

# Extract names from first columns
species_names <- species[[1]]
individual_names <- indivs[[1]]

# 6. Load phylogenetic trees
clade_tree <- read.tree(NUDNA_PHYLO)
gtrees_ind <- read.tree(GENETREES_NUDNA)

# 7. Create quartet table
b <- quartetTable(gtrees_ind, taxonnames = individual_names, random = 0)
b <- b[rowSums(b) != 4, colSums(b) != 0]

# 8. Create the quartet file for SNaQ
x <- b
y <- matrix(rep(0, dim(x)[1] * 8), ncol = 8, nrow = dim(x)[1])
colnames(y) <- c("t1", "t2", "t3", "t4", "CF12_34", "CF13_24", "CF14_23", "ngenes")

for(i in 1:dim(x)[1]) {
  count = 0
  for(j in 1:dim(x)[2]) {
    if(x[i,j] == 1) {
      count = count + 1
      y[i, count] = colnames(x)[j]
    }
  }
  y[i, "ngenes"] = x[i, dim(x)[2] - 2] + x[i, dim(x)[2] - 1] + x[i, dim(x)[2]]
  y[i, "CF12_34"] = as.numeric(x[i, dim(x)[2] - 2]) / as.numeric(y[i, "ngenes"])
  y[i, "CF13_24"] = as.numeric(x[i, dim(x)[2] - 1]) / as.numeric(y[i, "ngenes"])
  y[i, "CF14_23"] = as.numeric(x[i, dim(x)[2]]) / as.numeric(y[i, "ngenes"])
}

# 9. Write output files to the group-specific directory
setwd(group_dir)
write.csv(y, paste0(group, "_individualsCFs.csv"), row.names = FALSE)

# 10. Process the species tree
tree_new <- keep.tip(clade_tree, species_names)
tree_new$edge.length <- NULL
tree_new$node.label <- NULL

# 11. Write the final tree file
write.tree(tree_new, paste0(group, ".tre"))

# Print completion message
cat("Analysis completed for group:", group, "\n")
cat("Results saved to:", group_dir, "\n")
