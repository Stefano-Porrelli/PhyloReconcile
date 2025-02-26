#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#                   Script_01b - Concordance Analysis 1                        #
#                   SuperTRI analysis of phylogenetic signal                   #
#------------------------------------------------------------------------------#

# Set paths
BASE_DIR="$(pwd)/Phylogenomic_pipeline"
CONCORDANCE_DIR="${BASE_DIR}/03_Concordance_analyses"
SUPERTRI_OUT="${CONCORDANCE_DIR}/SuperTRI"
SUPERTRI_IN="${BASE_DIR}/01_initial_data/Input_files/SuperTRI_input"

# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate Phylogenomic_pipeline
# Check the Python version (>3 is needed for SuperTRI_v60)
PYTHON_VERSION=$(python -c 'import sys; print(sys.version_info[0])')
if [[ "$PYTHON_VERSION" -ne 3 ]]; then
    echo "Python version is not 3. Installing Python 3..."
    conda install -y python=3
else
    echo "Python 3 is already installed."
fi

# Ensure the required directories exist
mkdir -p "$CONCORDANCE_DIR"    
mkdir -p "$SUPERTRI_OUT"       

# Prepare input files for SuperTRI
cd "$SUPERTRI_OUT" || exit   # Change to SUPERTRI_OUT directory

# Retrieve SuperTRI v.60 (https://bitbucket.org/blaiseli/supertri/src/master/)
curl -O https://bitbucket.org/blaiseli/supertri/raw/master/supertri.py
chmod +x supertri.py

# Retrieve PAUP (MacOS)
curl -O http://phylosolutions.com/paup-test/paup4a168_osx.gz
gunzip paup4a168_osx.gz
chmod +x paup4a168_osx

# For Ubuntu/Linux, comment the MacOS block above and uncomment the block below
# curl -O http://phylosolutions.com/paup-test/paup4a169_ubuntu64.gz
# gunzip paup4a169_ubuntu64.gz
# chmod +x paup4a169_ubuntu64.gz

# Retrieve input data
cp "${SUPERTRI_IN}"/*.parts "${SUPERTRI_OUT}"   # Files describing bipartitions (MrBayes Output)
cp "${SUPERTRI_IN}"/*.tstat "${SUPERTRI_OUT}"   # Files of summary statistics (MrBayes Output)
cp "${SUPERTRI_IN}"/*.abs "${SUPERTRI_OUT}"    # Files with names of missing taxa for each marker (if applicable)
cp "${SUPERTRI_IN}/taxa.txt" "${SUPERTRI_OUT}"  # List of taxon names, in the same order as in the matrix used
cp "${SUPERTRI_IN}/genes.txt" "${SUPERTRI_OUT}" # List of marker names (prefixes)
cp "${SUPERTRI_IN}/paup_block.txt" "${SUPERTRI_OUT}" # PAUP block 

# Run SuperTRI 
python3 supertri.py --taxlist taxa.txt --datasets genes.txt --suffix .parts -o MRP.nex --root Panthera_leo_Ple1

# Run PAUP and generate synthesis tree with SBP support
./paup4a168_osx
exe paup_block.txt

# Map SuperTRI indices onto synthesis tree
# Will generate a file with 2 trees mapped with MPP and Nrep
python3 supertri.py --taxlist taxa.txt --datasets genes.txt --suffix .parts -o "${SUPERTRI_OUT}/MRP.nex" --root Panthera_leo_Ple1 --intree "${SUPERTRI_OUT}/synthesistree_boot.tre"

# Extract tree with Nrep metrics from SuperTRI analyses
awk '/tree repro/,/end;/' "${SUPERTRI_OUT}/synthesistree_boot.tre.withindices" | \
    sed -e 's/tree repro =//g' -e 's/end;//g' > "${SUPERTRI_OUT}/Nreps.tree"

# Extract tree with MPP metrics from SuperTRI analyses
awk '/tree meansup/,/end;/' "${SUPERTRI_OUT}/synthesistree_boot.tre.withindices" | \
    sed -e 's/tree meansup =//g' -e 's/end;//g' > "${SUPERTRI_OUT}/MPP.tree"

# Root trees
nw_reroot Nreps.tree Panthera_leo_Ple1 > Nreps_rooted.tree
nw_reroot MPP.tree Panthera_leo_Ple1 > MPP_rooted.tree

# Generate gCF/sCF for SuperTRI tree (which might differ from ML species tree)
iqtree2 -t Nreps_rooted.tree --gcf "${DATA_DIR}/ml_best.trees" -s "${CONCAT_ALIGNMENT}" --scf 100 --prefix 30AX_Nreps_tree --cf-verbose --df-tree --cf-quartet
iqtree2 -t MPP_rooted.tree --gcf "${DATA_DIR}/ml_best.trees" -s "${CONCAT_ALIGNMENT}" --scf 100 --prefix 30AX_MPP_tree --cf-verbose --df-tree --cf-quartet

# Cleanup
mkdir ./Result_trees
mv ./*.tree ./Result_trees
mv ./*.tre ./Result_trees
