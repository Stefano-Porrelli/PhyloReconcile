#!/bin/bash
#------------------------------------------------------------------------------#
#                             PhyloReconcile                                   #
#                                                                              #
#             Script_01d - Analysis of introgression with Dsuite               #
#                         (D-statistics, F-branch)                             #
#------------------------------------------------------------------------------#
# Set paths
BASE_DIR="$(pwd)/PhyloReconcile"
DATA_DIR="${BASE_DIR}/01_initial_data/30AX_MSAs"
CLADES_TREES="${BASE_DIR}/01_initial_data/Input_files/Dsuite_input/Clades_Trees"
CLADES_MSA="${BASE_DIR}/01_initial_data/Input_files/Dsuite_input/Clades_MSAs"
TAXASETS="${BASE_DIR}/01_initial_data/Input_files/Dsuite_input"
DSUITE_OUT="${BASE_DIR}/05_Dsuite_output"
DSUITE_DIR="$(pwd)/Dsuite"

# Ensure Conda is initialized properly
source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate the Conda environment
conda activate PhyloReconcile

# If not already installed, install Dsuite
if command -v "${DSUITE_DIR}/Build/Dsuite" &> /dev/null; then
    echo "Dsuite is already installed. Skipping installation."
else
    echo "Installing Dsuite..."
    git clone https://github.com/millanek/Dsuite.git "${DSUITE_DIR}"
    cd "${DSUITE_DIR}"
    make
    cd "${DSUITE_DIR}/utils"
    python3 setup.py install --user --prefix=
    cd "${BASE_DIR}"
fi

# Loop through each taxon set file
for TAXASET_FILE in "${TAXASETS}"/*_TaxaSet.txt; do
    # Extract prefix (e.g., "01_CETACEA")
    PREFIX=$(basename "${TAXASET_FILE}" _TaxaSet.txt)
    # Extract clade name (e.g., "CETACEA")
    CLADE_NAME=${PREFIX#*_}
    PREFIX_NUM=${PREFIX%_*}  # Get number prefix (e.g., "01")
    
    echo "Processing ${CLADE_NAME}..."
    
    # Create output directories
    mkdir -p "${DSUITE_OUT}"
    cd "${DSUITE_OUT}" || exit
    mkdir -p ./${PREFIX}/snp_output
    
    # Copy input files
    cp "${CLADES_MSA}/${PREFIX_NUM}_30AX_${CLADE_NAME}.fasta" ./${PREFIX}/snp_output
    cp "${DATA_DIR}/30AX_concatenated/30AX_partitions.txt" ./${PREFIX}/snp_output
    
    # Process the MSA files
    cd ./${PREFIX}/snp_output
    
    # Replace '?'s by 'N'
    sed -i'.bak' -e 's/\?/N/g' ${PREFIX_NUM}_30AX_${CLADE_NAME}.fasta
    # Replace 'N's by '-'
    sed -i'.bak' -e '/^>/!s/N/-/g' ${PREFIX_NUM}_30AX_${CLADE_NAME}.fasta
    # Remove temp files
    rm *.bak
    
    # Generate VCF with snp-sites
    snp-sites -v ${PREFIX_NUM}_30AX_${CLADE_NAME}.fasta > ${PREFIX_NUM}_30AX_${CLADE_NAME}.vcf
    
    # VCF parsing
    grep ^# ${PREFIX_NUM}_30AX_${CLADE_NAME}.vcf > header.txt
    grep -v "^#" ${PREFIX_NUM}_30AX_${CLADE_NAME}.vcf > snps.txt
    cut -d'=' -f2 30AX_partitions.txt | tr '-' '\t' > bed.txt
    
    # Filter VCF
    gzip ${PREFIX_NUM}_30AX_${CLADE_NAME}.vcf
    bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.2' -m2 -M2 -O z -o ${PREFIX_NUM}_30AX_${CLADE_NAME}_bcf.vcf.gz ${PREFIX_NUM}_30AX_${CLADE_NAME}.vcf.gz
    gzip -d -c ${PREFIX_NUM}_30AX_${CLADE_NAME}_bcf.vcf.gz > ${PREFIX_NUM}_30AX_${CLADE_NAME}_bcf.vcf
    grep ^# ${PREFIX_NUM}_30AX_${CLADE_NAME}_bcf.vcf > headerF.txt
    grep -v "^#" ${PREFIX_NUM}_30AX_${CLADE_NAME}_bcf.vcf > snpsF.txt
    rm ${PREFIX_NUM}_30AX_${CLADE_NAME}_bcf.vcf
    mv ${PREFIX_NUM}_30AX_${CLADE_NAME}.vcf.gz 30AXFULL_${CLADE_NAME}.vcf.gz
    mv ${PREFIX_NUM}_30AX_${CLADE_NAME}_bcf.vcf.gz 30AXFILTERED_${CLADE_NAME}.vcf.gz
    
    # Thinning (i.e., one SNP in 100bp window)
    vcftools --gzvcf 30AXFILTERED_${CLADE_NAME}.vcf.gz --thin 100 --recode --out 30AX_${CLADE_NAME}_thinned
    
    # Count thinned SNPs
    grep -v "^#" 30AX_${CLADE_NAME}_thinned.recode.vcf > snps_thin.txt
    mv 30AX_${CLADE_NAME}_thinned.recode.vcf 30AX_${CLADE_NAME}_thinned.vcf
    gzip 30AX_${CLADE_NAME}_thinned.vcf
    
    # SNP statistics
    initial=$(wc -l < snps.txt)
    filtered=$(wc -l < snpsF.txt)
    thinned=$(wc -l < snps_thin.txt)
    echo -e "Initial\t${initial}\nFiltered\t${filtered}\nThinned\t${thinned}" > SNPstat_${CLADE_NAME}.txt
    
    # Take the filtered vcf for downstream analysis
    cp 30AXFILTERED_${CLADE_NAME}.vcf.gz "${DSUITE_OUT}/${PREFIX}/30AXDsuite.vcf.gz"
    # Take the thinned vcf for downstream analysis (faster)
    # cp 30AX_${CLADE_NAME}_thinned.vcf.gz "${DSUITE_OUT}/${PREFIX}/30AXDsuite.vcf.gz"

    # Prepare files for Dsuite/Dtrios
    cd "${DSUITE_OUT}" || exit
    
    # Copy the tree file
    cp "${CLADES_TREES}/TREE_${CLADE_NAME}.tree" "${DSUITE_OUT}/${PREFIX}"
    
    # Copy the taxa list file
    cp "${TAXASETS}/${PREFIX}_TaxaSet.txt" "${DSUITE_OUT}/${PREFIX}"
    
    cd ./${PREFIX}
    
    # Run Dsuite/Dtrios
    "${DSUITE_DIR}/utils/DtriosParallel" -t TREE_${CLADE_NAME}.tree ${PREFIX}_TaxaSet.txt 30AXDsuite.vcf.gz
    
    # Plot F-branch using species tree as guide (p=0.05)
    "${DSUITE_DIR}/Build/Dsuite" Fbranch -p 0.05 TREE_${CLADE_NAME}.tree DTparallel_${PREFIX}_TaxaSet_combined_tree.txt > fbranch_p_005.txt
    python3 "${DSUITE_DIR}/utils/dtools.py" -n fbranch_p_005 fbranch_p_005.txt TREE_${CLADE_NAME}.tree
    
    # Plot F-branch using species tree as guide
    "${DSUITE_DIR}/Build/Dsuite" Fbranch TREE_${CLADE_NAME}.tree DTparallel_${PREFIX}_TaxaSet_combined_tree.txt > fbranch.txt
    python3 "${DSUITE_DIR}/utils/dtools.py" fbranch.txt TREE_${CLADE_NAME}.tree
    
    echo "Completed analysis for ${CLADE_NAME}"
    echo "----------------------------------------"
done

echo "All analyses completed!"