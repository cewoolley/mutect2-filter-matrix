#!/bin/bash

# Call filter-VCF.py to create the matrix for all samples using 185 drivers from 100KGP.

# Directory containing VCF files
VCF_DIR="/exports/igmm/eddie/tomlinson-Polyp-WGS-RNA/GRAMPIAN/grampian_VEP/tumour_normal"
# Output TSV file
OUTPUT_MATRIX="/exports/igmm/eddie/tomlinson-Polyp-WGS-RNA/GRAMPIAN/grampian_VEP/tumour_normal/grampian_tumour_normal_185_driver_matrix.tsv"
GENE_LIST="/exports/igmm/eddie/tomlinson-Polyp-WGS-RNA/GRAMPIAN/grampian_VEP/tumour_normal/185drivers.txt"

# Create or clear the output file
echo -e "Sample\t$(cat $GENE_LIST | tr '\n' '\t')" > "$OUTPUT_MATRIX"

# Iterate over each VCF file in the directory
for vcf_file in "$VCF_DIR"/*.vep.vcf; do
    echo "Processing $vcf_file..."
    # Call the Python script to create the mutation matrix
    python filter-VCF.py --matrix_only --gene_list "$GENE_LIST" --matrix_output temp_matrix.tsv "$vcf_file"
    
    # Append the temp matrix to the output matrix, skipping the header
    tail -n +2 temp_matrix.tsv >> "$OUTPUT_MATRIX"
    
    # Clean up the temporary matrix file
    rm temp_matrix.tsv
done

echo "Combined mutation matrix saved to $OUTPUT_MATRIX"