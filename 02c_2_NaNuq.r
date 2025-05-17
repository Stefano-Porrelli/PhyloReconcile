#!/usr/bin/env Rscript
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#                    Script_02c - 2. Network analysis SNaQ                     #
#                                                                              # 
#------------------------------------------------------------------------------#

# 1. Load required libraries
library('MSCquartets')
library('ape')

# 2. Set parameters
group = " " # Name of group and Directory
alpha = 0.05 # Significance/alpha level 

# 3. Setup directory structure
BASE_DIR <- file.path(getwd(), "PhyloReconcile")
NUDNA_PHYLO <- file.path(BASE_DIR, "02_phylogenies", "01b_30AX_ML_ASTRAL", "30AX_ASTRAL_ML_species.tree")
GENETREES_NUDNA <- file.path(BASE_DIR, "02_phylogenies", "01b_30AX_ML_ASTRAL", "ml_best.trees")  # Using NUDNA as specified

MTDNA_PHYLO <- file.path(BASE_DIR, "02_phylogenies", "02a_mtDNA_phylogeny", "trees", "MtDNA_concat.tree" )
GENETREES_MTDNA <- file.path(BASE_DIR, "02_phylogenies", "02a_mtDNA_phylogeny", "MtDNA_loci.treefile")  # Using NUDNA as specified


INPUTS <- file.path(BASE_DIR, "01_initial_data", "input_files", "SNaQ_NANUQ")
MSC_RESULTS <- file.path(BASE_DIR, "09_Network_Analyses")

# 4. Create/navigate to results directory based on group name
group_dir <- file.path(MSC_RESULTS, group)
if (!dir.exists(group_dir)) {
  dir.create(group_dir, recursive = TRUE)
}

# 5. Load input files
setwd(INPUTS)
species <- read.delim(paste0(group, "_species.txt"), header = FALSE)
outgroup <- species[length(species)]

# Extract names from first column
species_names <- species[[1]]

# 6. Load gene trees
genedata_ind = GENETREES_NUDNA  # Using NUDNA gene trees
gtrees_ind <- read.tree(genedata_ind)

# 7. Create quartet table for NANUQ
b <- quartetTable(gtrees_ind, taxonnames = species_names, random = 0)
b <- b[rowSums(b) != 4, colSums(b) != 0]  # Edit it to look the way we want for SNaQ

# 8. Move to output directory and write quartet file
setwd(group_dir)
write.csv(b, paste0(group, "_fullspeciesCFs.csv"), row.names = FALSE)

# 9. Read and modify the file for NANUQ format
u <- read.csv(paste0(group, "_fullspeciesCFs.csv"))
colnames(u)[dim(u)[2] - 2] = "12|34"
colnames(u)[dim(u)[2] - 1] = "13|24"
colnames(u)[dim(u)[2]] = "14|23"

# 10. Run NANUQ analysis
z <- NANUQ(as.matrix(u), 
           outfile = paste0(group), 
           alpha = alpha,  # Using the variable rather than string "alpha"
           beta = 0.95, 
           plot = TRUE)

# Print completion message
cat("NANUQ analysis completed for group:", group, "\n")
cat("Results saved to:", group_dir, "\n")
