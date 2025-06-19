# PhyloReconcile <!-- omit in toc --> 
A phylogenomic pipeline for resolving species relationships through multilocus data, accounting for discordance, incomplete lineage sorting (ILS), and introgression.

<img src="01_initial_data/Pipeiline_flowchart/PhyloReconcile_flowchart.pdf" alt="Pipeline Workflow" width="600">

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
  - [4.7 Script `01d_Dsuite_analyses.sh`](#47-script-01d_dsuite_analysessh)
  - [4.8 Script `02a_1_mtDNA_tree.sh`](#48-script-02a_1_mtdna_treesh)
  - [4.9 Script `02a_2_MSY_tree.sh`](#49-script-02a_2_msy_treesh)
  - [4.10 Script `02a_3_Mitonuclear_tanglegrams.sh`](#410-script-02a_3_mitonuclear_tanglegramssh)
  - [4.11 Script `02b_1_Phyparts_mtDNA.sh`](#411-script-02b_1_phyparts_mtdnash)
  - [4.12 Script `02b_2_DiscoVista_mtDNA.sh`](#412-script-02b_2_discovista_mtdnash)
  - [4.13 Script `02c_1_SNaQ_prep.r`](#413-script-02c_1_snaq_prepr)
  - [4.14 Script `02c_2_NaNuq.r`](#414-script-02c_2_nanuqr)
  - [4.15 Script `02c_3_SNaQ_Networks.jl`](#415-script-02c_3_snaq_networksjl)

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
├── 06_comparison_inheritance_modes/
├── 07_Phyparts_MtDNA_MSY/
├── 08_DiscoVista_mtDNA/
└── 09_Network_Analyses/
```
This structure is designed to organize the various analyses and outputs in a logical manner:

- 01_initial_data: Contains the input multiple sequence alignments and other required input files
- 02_phylogenies: Stores the gene trees and species trees from different methods
- 03_Concordance_analyses: Contains concordance analyses results
- 04_DiscoVista_analysis: Stores the outputs from DiscoVista discordance visualization
- 05_Dsuite_output: Contains introgression analysis results from Dsuite
- 06_comparison_inheritance_modes: Stores tanglegram comparisons between different inheritance modes
- 07_Phyparts_MtDNA_MSY: Contains Phyparts analyses for mitochondrial and Y-chromosome data
- 08_DiscoVista_mtDNA: Stores DiscoVista analyses for mitochondrial data
- 09_Network_Analyses: Contains phylogenetic network analyses (SNaQ, NANUQ)

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

**For SNaQ/NANUQ network analyses, additional files are required:**

- Individual-to-species mapping files
- Species lists for each clade
- Individual lists for each clade
  
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
- `02a_3_Mitonuclear_tanglegrams.sh`: Assess and visualize topological discordances between phylogenies generated from markers of different inheritance mode (biparental, maternal and paternal)
- `02b_1_Phyparts_mtDNA.sh`: Performs Phyparts analyses for mitochondrial and Y-chromosome data
- `02b_2_DiscoVista_mtDNA.sh`: Performs DiscoVista analyses for mitochondrial data
- `02c_1_SNaQ_prep.r`: Prepares data for phylogenetic network analysis
- `02c_2_NaNuq.r`: Performs NANUQ network analysis
- `02c_3_SNaQ_Networks.jl`: Performs SNaQ phylogenetic network inference

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

These visualizations can be found in `04_DiscoVista_analysis/results/` directory, organized by taxonomic level and specific nodes of interest.

### 4.7 Script `01d_Dsuite_analyses.sh`

**Analysis of introgression with Dsuite**

Execute the script to perform D-statistics and F-branch analyses for introgression detection:

`./01d_Dsuite_analyses.sh`

The script:

- Installs Dsuite if not already present
- Processes multiple sequence alignments to generate VCF files
- Filters SNPs to remove sites with high missing data
- Performs SNP thinning (one SNP per 100bp window)
- Runs Dsuite Dtrios analysis for each clade
- Generates F-branch plots using species trees as guides
- Creates statistical summaries and visualizations

The analyses are performed separately for each taxonomic group defined in the taxon set files. Results include D-statistics tables, F-branch plots, and significance tests for introgression.

The outputs are stored in `05_Dsuite_output/` with separate directories for each analyzed clade.

### 4.8 Script `02a_1_mtDNA_tree.sh`

**Generation of Mitochondrial Phylogeny**

Execute the script to generate mitochondrial phylogenies using concatenation and MSC approaches:

`./02a_1_mtDNA_tree.sh`

The script:

- Infers ML mitochondrial tree using concatenation method with IQ-TREE2
- Generates individual mitochondrial gene/locus trees
- Computes gCF and sCF for the concatenated ML mitochondrial tree
- Runs ASTRAL with multilocus bootstrap for MSC analysis
- Runs ASTRAL with local posterior probabilities
- Computes gCF and sCF for both ASTRAL trees

The outputs are stored in `02_phylogenies/02a_mtDNA_phylogeny/trees/`

### 4.9 Script `02a_2_MSY_tree.sh`

**Generation of MSY Phylogeny (AMELY/DDX3Y)**

Execute the script to validate Y chromosome sequences and generate Y chromosome phylogenies:

`./02a_2_MSY_tree.sh`

The script:

- Compares X and Y chromosome sequences for AMEL and DDX3 genes to validate Y chromosome sequence authenticity
- Generates concatenated reference tree for MSY dataset
- Generates individual Y chromosome loci trees
- Computes gCF and sCF for the concatenated ML MSY tree
- Runs ASTRAL analyses (multilocus bootstrap and local PP)
- Computes gCF and sCF for ASTRAL MSY trees

The X/Y comparison trees help validate the authenticity of Y chromosome sequences by ensuring they cluster separately from their X homologs.

The outputs are stored in `02_phylogenies/02b_MSY_phylogeny/trees/`

### 4.10 Script `02a_3_Mitonuclear_tanglegrams.sh`

**Comparison of topological discordance between inheritance modes**

Execute the script to assess and visualize topological discordances between phylogenies from different inheritance modes:

`./02a_3_Mitonuclear_tanglegrams.sh`

The script generates four R scripts that create tanglegram comparisons:

- **AMELY vs DDX3Y comparison**: Compares the two Y chromosome markers
- **Mitochondrial vs Nuclear DNA comparison**: Compares maternal vs biparental inheritance
- **Mitochondrial vs MSY DNA comparison**: Compares maternal vs paternal inheritance  
- **MSY vs Nuclear DNA comparison**: Compares paternal vs biparental inheritance

Each tanglegram shows:
- Side-by-side tree topologies with connecting lines between corresponding taxa
- Node support values below 100 are displayed
- Curved connecting lines colored to show topological agreements/disagreements

The outputs are PDF files stored in `06_comparison_inheritance_modes/`

### 4.11 Script `02b_1_Phyparts_mtDNA.sh`

**Detection of ILS with PhyParts - Mitochondrial and Y chromosome data**

Execute the script to perform Phyparts analyses on organellar genomes:

`./02b_1_Phyparts_mtDNA.sh`

The script:

- Creates a Python 2.7 environment for Phyparts compatibility
- Installs Phyparts using Maven
- Processes mitochondrial gene trees and species tree
- Roots all trees using the specified outgroup
- Runs Phyparts analysis for mitochondrial data (14 gene trees)
- Processes Y chromosome gene trees and species tree  
- Runs Phyparts analysis for MSY data (2 gene trees)
- Generates pie chart visualizations for both analyses

The pie charts show the proportion of gene trees that support, conflict with, or are uninformative for each branch in the species tree.

The outputs are stored in `07_Phyparts_MtDNA_MSY/`

### 4.12 Script `02b_2_DiscoVista_mtDNA.sh`

**Visualization of gene tree and species tree discordance with DiscoVista - Mitochondrial data**

Execute the script to perform DiscoVista analyses specifically for mitochondrial data:

`./02b_2_DiscoVista_mtDNA.sh`

The script:

- Processes mitochondrial gene trees into individual directories for DiscoVista
- Organizes species trees for comparative analysis
- Generates clade definition files for mitochondrial-specific analyses
- Runs DiscoVista discordance analysis on gene trees showing proportion of gene trees supporting/rejecting specific clades
- Performs relative frequency analysis to determine topology frequencies around focal branches
- Uses Docker to run DiscoVista analyses

The analyses focus on clade monophyly patterns in mitochondrial data and can reveal organelle-specific phylogenetic patterns.

The outputs are stored in `08_DiscoVista_mtDNA/results/`

### 4.13 Script `02c_1_SNaQ_prep.r`

**Data preparation for SNaQ phylogenetic network analysis**

Execute the R script to prepare quartet concordance factor data for network analysis:

`Rscript 02c_1_SNaQ_prep.r`

The script:

- Loads required R libraries (MSCquartets, ape, gtools, stringr)
- Reads individual and species mapping files for the specified group
- Loads nuclear DNA gene trees and species tree
- Creates quartet table from gene trees using MSCquartets
- Formats data for SNaQ input with concordance factors (CF12_34, CF13_24, CF14_23)
- Processes species tree to include only taxa of interest
- Writes output files for downstream network analysis

**Required input files**:
- `[group]_individuals.txt`: List of individuals in the clade
- `[group]_species.txt`: List of species in the clade  
- `[group]_map.csv`: Mapping of individuals to species

**Note**: Edit the script to specify the group name before running.

The outputs are stored in `09_Network_Analyses/[group]/`

### 4.14 Script `02c_2_NaNuq.r`

**Network analysis with NANUQ**

Execute the R script to perform NANUQ network analysis:

`Rscript 02c_2_NaNuq.r`

The script:

- Loads MSCquartets and ape libraries
- Reads species lists and gene trees for the specified group
- Creates quartet table for NANUQ analysis
- Runs NANUQ analysis to test for network-like evolution
- Performs statistical tests with specified alpha level (default 0.05)
- Generates plots and statistical summaries

NANUQ tests whether quartet concordance factors are consistent with a tree model or suggest network-like evolution patterns.

**Note**: Edit the script to specify the group name and alpha level before running.

The outputs are stored in `09_Network_Analyses/[group]/`

### 4.15 Script `02c_3_SNaQ_Networks.jl`

**Phylogenetic Network Analysis Using PhyloNetworks in Julia**

Execute the Julia script to perform SNaQ phylogenetic network inference:

`julia 02c_3_SNaQ_Networks.jl`

The script:

- Loads PhyloNetworks, PhyloPlots, CSV, DataFrames, RCall, and QuartetNetworkGoodnessFit
- Maps individual-level concordance factors to species-level
- Runs SNaQ analysis with increasing hybridization levels (h=0 to h=5)
- Performs model selection based on pseudo-likelihood scores
- Calculates goodness-of-fit statistics for each network
- Generates SVG visualizations of inferred networks
- Tests statistical significance of reticulation events
- Creates results summary with network topologies and statistics

**Prerequisites**: 
- Julia installation with PhyloNetworks.jl package
- R integration through RCall

**Note**: Edit the script to specify the group name and outgroup before running.

The script systematically searches for the optimal number of hybridization events and provides statistical support for inferred reticulation patterns.

The outputs include network files, visualizations, and statistical summaries stored in `09_Network_Analyses/[group]/`
