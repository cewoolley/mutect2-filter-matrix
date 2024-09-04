#!/bin/bash
for vcf in *filtered.vcf.gz
do
	echo "Filtering $vcf for PASS flags and min AD of 5 in >=1 tumour sample. See new header lines for info."
	bn=$(basename $vcf .vcf.gz)
	python filter-VCF.py "$vcf" --output_file "${bn}.PASS.vcf.gz" --min_alt_ad=5
done
