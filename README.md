# PhyloReconcile <!-- omit in toc --> 

A phylogenomic pipeline for resolving species relationships through multilocus data, accounting for discordance, incomplete lineage sorting (ILS), and introgression.

## Table of contents <!-- omit in toc --> 
- [1. Introduction](#1-introduction)
- [2. Preparing software for analyses](#2-preparing-software-for-analyses)
  - [2.1. Installation of required software](#21-installation-of-required-software)
  - [2.2. Directory structure](#22-directory-structure)
- [3. Preparing data for analyses](#3-preparing-data-for-analyses)
  - [3.1 Required input files](#31-required-input-files)
  - [3.2. Input MSAs](#32-input-msas)
- [4. Running the pipeline](#4-running-the-pipeline)
  - [4.1 Script `00_PhyloReconcile_setup.sh`](#41-script-00_phyloreconcile_setupsh)
  - [4.2 Script `01a_ASTRAL_speciestree_MSC.sh`](#42-script-01a_astral_speciestree_mscsh)
  - [4.3 Script `01b_Concordance_analyses_1_SuperTRI.sh`](#43-script-01b_concordance_analyses_1_supertrish)
  - [4.4 Script `01b_Concordance_analyses_2_CF_DF.sh`](#44-script-01b_concordance_analyses_2_cf_dfsh)
  - [4.5 Script `01b_Concordance_analyses_3_Phyparts.sh`](#45-script-01b_concordance_analyses_3_phypartssh)
  - [4.6 Script `01c_DiscoVista_analyses.sh`](#46-script-01c_discovista_analysessh)

## 1. Introduction
PhyloReconcile is a phylogenomic pipeline for reconstructing species trees from multilocus datasets, addressing gene tree discordance, incomplete lineage sorting (ILS), and introgression through robust analyses with Maximum Likelihood (ML) and Bayesian Inference (BI) methods, multi-species coalescent (MSC) model, discordance visualization, and introgression detection. The pipeline consists of a set of BASH and R scrips for UNIX environment such as Linux or MacOS.

The pipeline integrates several phylogenetic tools including IQ-TREE2, ASTRAL, MrBayes, DiscoVista and Dsuite to provide a comprehensive analysis of phylogenetic relationships while accounting for possible sources of discordance.

## 2. Preparing software for analyses

Software required to run the PhyloReconcile can be installed by executing the script  `00_PhyloReconcile_setup.sh`, assuming [Miniconda](https://docs.anaconda.com/miniconda/install/) is already installed. If Miniconda is not installed, lines 23-28 of the script `00_PhyloReconcile_setup.sh` can be uncommented (remove #). This will install Miniconda when executing the script.

Before running PhyloReconcile, the scripts need to be downloaded into `$HOME` directory and made executable by running the command `chmod +x  *.sh` from the `$HOME` directory.

### 2.1. Installation of required software

Executing the script `00_PhyloReconcile_setup.sh` will install all necessary software and create the main directory and subdirectories. The script creates a Conda environment named "PhyloReconcile" with the following tools:

- IQ-TREE2: For maximum likelihood (ML) phylogenetic inference;
- ASTRAL: For species tree inference under the multi-species coalescent (MSC) model;
- MrBayes: For Bayesian inference (BI) of phylogenies;
- R and packages: For statistical analyses and plotting (viridis, ggplot2, dplyr, ggrepel, ggally, entropy, patchwork);
- Additional tools: snp-sites, samtools, vcftools, bcftools, and newick_utils.

Executing the script `00_PhyloReconcile_setup.sh` will also retrieve example data such as MSAs and other input files required by some of the software (such as Dsuite and DiscoVista).

To run the script, make it executable by running the command `chmod +x  *.sh` from the `$HOME` directory (if not done already) and run the script with `./00_PhyloReconcile_setup.sh`

### 2.2. Directory structure
The script creates the following directory structure:

```plaintext
PhyloReconcile/
PhyloReconcile/
├── 01_initial_data/
│   ├── 30AX_MSAs/
│   ├── mtDNA_MSA/
│   ├── YDNA_MSA/
│   └── Input_files/
│       ├── SuperTRI_input/
│       ├── Dsuite_input/
│       └── DiscoVista_parameters/
├── 02_phylogenies/
│   ├── 01a_30AX_MLgenetrees/
│   ├── 01b_30AX_ML_ASTRAL/
│   ├── 01c_30AX_BIgenetrees/
│   ├── 01d_30AX_BI_ASTRAL/
│   ├── 01e_30AX_MLconcat/
│   ├── 02a_mtDNA_phylogeny/
│   └── 02b_MSY_phylogeny/
├── 03_Concordance_analyses/
│   ├── Concordance_factors/
│   ├── SuperTRI/
│   ├── CF_DF_analysis/
│   └── Phyparts_quartet_analyses/
├── 04_DiscoVista_analysis/
├── 05_Dsuite_output/
└── 06_comparison_inheritance_modes/ 
└── ...
```
This structure is designed to organize the various analyses and outputs in a logical manner:

- 01_initial_data: Contains the input multiple sequence alignments and other required input files
- 02_phylogenies: Stores the gene trees and species trees from different methods
- 03_Concordance_analyses: Contains concordance analyses results
- 05_DiscoVista_analysis: Stores the outputs from DiscoVista discordance visualization

## 3. Preparing data for analyses

### 3.1 Required input files

Before running the pipeline, you need to prepare the following input files:

**Autosomal data**:

- Multiple Sequence Alignments (MSAs): For each gene or locus in NEXUS format
- Concatenated alignment: All loci concatenated into a single alignment
- Partition file: Defining the boundaries of each locus in the concatenated alignment
- MrBayes blocks: NEXUS files with MrBayes command blocks for each gene

**Mitochondrial data**:

- MSA for mitochondrial genes in NEXUS format
- Partition file for mitochondrial genes

**Y-chromosome data**:

- MSAs for Y-chromosome markers (AMELY, DDX3Y) in NEXUS format
- MSA for Y-chromosome genes and their X-homologs in NEXUS format
- Partition file for Y-chromosome genes

**For concordance analyses, additional files are required**:

- SuperTRI input files: Files describing bipartitions (.parts), summary statistics (.tstat), taxa lists (taxa.txt), gene lists (genes.txt), and PAUP blocks

**For DiscoVista analyses, additional annotation files are required:**

- Taxonomic annotations at different levels (family/subfamily, tribe/subtribe)
- Clade definition files
- Model order files

**For Dsuite analyses, additional files are required:**

- Taxa set files: Lists of taxa to include in the analysis
- Tree files: Guide trees for each clade of interest
  
### 3.2. Input MSAs

The input Multiple Sequence Alignments (MSAs) should be placed in the appropriate directories:

- Autosomal MSAs: `01_initial_data/30AX_MSAs/`
- Mitochondrial MSAs: `01_initial_data/mtDNA_MSAs/`
- Y-chromosome MSAs: `01_initial_data/YDNA_MSAs/`

Each MSA should represent a single gene or locus and be in NEXUS format.
For autosomal data, a concatenated alignment of all loci should be placed in `01_initial_data/30AX_MSAs/30AX_concatenated/30AX_concatenated.fasta` with a corresponding partition file at `01_initial_data/30AX_MSAs/30AX_concatenated/30AX_partitions.txt`.

## 4. Running the pipeline

The PhyloReconcile pipeline consists of a series of scripts that should be executed in order:

- `00_PhyloReconcile_setup.sh`: Sets up the environment and installs required software
- `01a_ASTRAL_speciestree_MSC.sh`: Generates gene trees and species trees from autosomal nuclear markers
- `01b_Concordance_analyses_1_SuperTRI.sh`: Performs SuperTRI analysis of phylogenetic signal
- `01b_Concordance_analyses_2_CF_DF.sh`: Performs analysis of concordance and discordance factors
- `01b_Concordance_analyses_3_Phyparts.sh`: Performs quartet-based analyses and visualizes quartet discordance with pie charts in Phyparts
- `01c_DiscoVista_analyses.sh`: Visualizes gene tree and species tree discordance
- `01d_Dsuite_analyses.sh`: Performs analysis of introgression with Dsuite
- `02a_1_mtDNA_tree.sh`: Generates mitochondrial phylogeny
- `02a_2_MSY_tree.sh`: Validates authenticity of Y chromosome sequences and generates Y chromosome phylogeny
- `02a_3_Mitonucleat_tanglegrams.sh`: Assess and visualize topological discordances between phylogenies generated from markers of different inheritance mode (biparental, maternal and paternal)

### 4.1 Script `00_PhyloReconcile_setup.sh`
This script sets up the environment and software required for the pipeline. It creates the necessary directory structure, installs required software through Conda, and prepares the environment for subsequent analyses.

Execute with:

`./00_PhyloReconcile_setup.sh`

This will:

- Create the directory structure
- Install all required software via Conda
- Install additional tools (snp-sites, samtools, vcftools, bcftools, newick_utils)

### 4.2 Script `01a_ASTRAL_speciestree_MSC.sh`

**Generating gene trees and species trees**

Execute `01a_ASTRAL_speciestree_MSC.sh` to generate gene trees using both Maximum Likelihood (IQ-TREE2) and Bayesian Inference (MrBayes), and species trees using ASTRAL under the multi-species coalescent model:

`./01a_ASTRAL_speciestree_MSC.sh`

This script performs the following steps:

- **ML gene trees**: Infers Maximum Likelihood gene trees using IQ-TREE2 with 1000 ultrafast bootstrap replicates
- **ML species tree**: Uses ASTRAL to infer a species tree from the ML gene trees
- **Concordance factors**: Calculates gene concordance factors (gCF) and site concordance factors (sCF) for the ML species tree
- **BI gene trees**: Infers Bayesian gene trees using MrBayes
- **BI species tree**: Uses ASTRAL to infer a species tree from the BI gene trees
- **ML concatenated tree**: Infers a Maximum Likelihood tree from the concatenated alignment

The outputs are stored in their respective directories under `02_phylogenies/`

### 4.3 Script `01b_Concordance_analyses_1_SuperTRI.sh`

**SuperTRI analyses**

Execute the script to assess the level of agreement between gene trees and species trees with SuperTRI:

`./01b_Concordance_analyses_1_SuperTRI.sh`

The script:

- Retrieves SuperTRI v.60 and PAUP (for either MacOS or Linux)
- Prepares input files for SuperTRI analysis
- Runs SuperTRI to generate synthesis trees
- Processes the output to extract trees with different support metrics (Nreps, MPP, SBP)
- Calculates concordance factors (gCF/sCF) for the SuperTRI trees

Output trees are stored in `03_Concordance_analyses/SuperTRI/Result_trees/`

### 4.4 Script `01b_Concordance_analyses_2_CF_DF.sh`

Analyzes concordance factors (CF) and discordance factors (DF) to test for incomplete lineage sorting (ILS) and generates a series of plots visualizing the relationships between different metrics.
Execute with:

`./01b_Concordance_analyses_2_CF_DF.sh`

The script:

- Prepares data from the previous analyses
- Creates an R script to analyze the concordance and discordance factors
- Calculates correlations between various metrics (gCF, sCF, bootstrap support, branch length)
- Generates plots visualizing these relationships

The outputs, including plots and statistics, are stored in `03_Concordance_analyses/CF_DF_analysis/`

### 4.5 Script `01b_Concordance_analyses_3_Phyparts.sh`
This script performs phylogenetic analysis using Phyparts to visualize gene tree support for the species tree.
Execute with:

`./01b_Concordance_analyses_3_Phyparts.sh`

The script:

- Creates a Python 2.7 environment with the necessary dependencies
- Installs Phyparts (requires Maven)
- Prepares and roots gene trees
- Runs Phyparts to analyze gene tree support for the species tree
- Generates pie charts visualizing the support for each branch

The outputs are stored in `03_Concordance_analyses/Phyparts_quartet_analyses/`

### 4.6 Script `01c_DiscoVista_analyses.sh`
This script uses DiscoVista to visualize gene tree and species tree discordance at different taxonomic levels.
Execute with:

`./01c_DiscoVista_analyses.sh`

This script uses DiscoVista to generate several visualizations:

- **Species tree discordance:** Compares clade support between different species tree inference methods
- **Gene tree discordance:** Shows the proportion of gene trees that support or reject specific clades
- **Relative frequency analysis:** Determines the frequency of alternative topologies around focal branches

The analyses are performed at different taxonomic levels (family/subfamily, tribe/subtribe) and for specific nodes of interest in the phylogeny.

DiscoVista generates several visualization types:

- Clade support heatmaps: Show support for specific clades across different species tree methods
- Gene tree support barplots: Show the proportion of gene trees that support, reject, or are indecisive about specific clades
- Relative frequency plots: Show the frequency of alternative topologies around focal branches

These visualizations can be found in `05_DiscoVista_analysis/results/` directory, organized by taxonomic level and specific nodes of interest.
