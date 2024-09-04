# filter-VCF.py
# connor.woolley@oncology.ox.ac.uk 2024
# Further filter Mutect2 VCF files and/or creates a mutation matrix across a list of genes of interest.
#
# Mutect2 variant calling -> GATK filterMutectCalls -> filter-VCF.py (this script) -> mutation matrix (optional!)
#
#
# Key functionalities:
# 1. Filters VCF files based on specified criteria (e.g., minimum alternate allele depth --min_alt_ad)
# 2. Processes tumour samples or all samples as specified
# 3. Supports different filtering modes: 'any', 'all', or 'majority'
# 4. Creates a mutation matrix from filtered VCF files
# 5. Generates detailed logging of the filtering process per-sample.
#
# The script can be run in two modes!
# a) Filtering mode: Filters the input VCF and optionally creates a mutation matrix (if already VEP annotated!)
# b) Matrix-only mode: Creates a mutation matrix without filtering the VCF
#
# Usage examples:
# python filter-VCF.py input.vcf --output_file filtered.vcf --min_alt_ad 5
# python filter-VCF.py input.vep.vcf --matrix_only --gene_list genes.txt --matrix_output matrix.tsv
#
# Variants can be filtered based on the following criteria:
# --filter_mode 'any' (default)
# --filter_mode 'all'
# --filter_mode 'majority'
#
# Info on modes:
# 'any': Keep the variant if it passes the filter in any of the processed samples.
#        This is the most lenient mode, allowing variants to be retained if they meet
#        the criteria in at least one sample.
# 'all': Keep the variant only if it passes the filter in all processed samples.
#        This is the strictest mode, requiring the variant to meet the criteria
#        across all samples to be retained.
# 'majority': Keep the variant if it passes the filter in more than half of the processed samples.
#             This mode strikes a balance between 'any' and 'all', requiring the variant
#             to meet the criteria in most samples but allowing for some exceptions.



import argparse
import gzip
import re
from collections import Counter, defaultdict
import datetime
import os
import csv

def find_tumor_samples(header_lines):
    tumor_samples = []
    for line in header_lines:
        match = re.search(r'##(tumor|tumour)_sample=(\S+)', line) # We look for tumor or tumour. Mutect2 should spell as tumor, as will this script in variables for consistency!
        if match:
            tumor_samples.append(match.group(2))
    return tumor_samples

def get_filter_mode_description(filter_mode):
    if filter_mode == 'any':
        return "Variant kept if it passes in any processed sample"
    elif filter_mode == 'all':
        return "Variant kept if it passes in all processed samples"
    elif filter_mode == 'majority':
        return "Variant kept if it passes in more than half of processed samples"
    else:
        return "Unknown filter mode"

def filter_vcf(input_file, output_file, min_alt_ad, use_tumor_samples=True, specific_sample=None, filter_mode='any'):
    stats = {
        'total_variants': 0,
        'passed_variants': 0,
        'filtered_variants': 0,
        'sample_pass_counts': Counter(),
        'chrom_counts': Counter(),
        'chrom_sample_counts': defaultdict(lambda: defaultdict(int)),
        'filter_reasons': Counter(),
    }
    
    with (gzip.open(input_file, 'rt') if input_file.endswith('.gz') else open(input_file, 'r')) as infile, \
         (gzip.open(output_file, 'wt') if output_file.endswith('.gz') else open(output_file, 'w')) as outfile:
        
        # Initialise empty variables
        header_lines = []
        samples = []
        tumor_samples = []
        tumor_indices = []

        for line in infile:
            if line.startswith('##'):
                header_lines.append(line)
                outfile.write(line)
                continue
            
            if line.startswith('#CHROM'):
                # Add summary info to the header
                summary = f"##FILTER_SUMMARY=<Script=VCF_PASS_AD_filter,MinAltAD={min_alt_ad},FilterMode={filter_mode},Date={datetime.datetime.now().strftime('%Y-%m-%d')}>\n"
                outfile.write(summary)
                
                # Add filter mode explanation
                filter_explanation = f"##FILTER_MODE_EXPLANATION=<{get_filter_mode_description(filter_mode)}>\n"
                outfile.write(filter_explanation)
                
                header_fields = line.strip().split('\t')
                samples = header_fields[9:]
                
                if use_tumor_samples:
                    tumor_samples = find_tumor_samples(header_lines)
                    if tumor_samples:
                        tumor_indices = [samples.index(sample) for sample in tumor_samples if sample in samples]
                        if not tumor_indices:
                            raise ValueError(f"No tumor samples found in VCF samples. Tumor samples in header: {', '.join(tumor_samples)}, VCF samples: {', '.join(samples)}")
                    else:
                        raise ValueError("No tumor samples found in VCF header. Use --all_samples to process all samples.")
                elif specific_sample:
                    if specific_sample in samples:
                        tumor_indices = [samples.index(specific_sample)]
                    else:
                        raise ValueError(f"Specified sample '{specific_sample}' not found in VCF samples: {', '.join(samples)}")
                else:
                    tumor_indices = list(range(len(samples)))
                
                outfile.write(line)
                stats['processed_samples'] = [samples[i] for i in tumor_indices]
                for sample in stats['processed_samples']:
                    stats['sample_pass_counts'][sample] = 0
                continue
            
            stats['total_variants'] += 1
            fields = line.strip().split('\t')
            chrom, pos, id_, ref, alt = fields[:5]
            filter_field = fields[6]
            format_field = fields[8]
            sample_fields = fields[9:]
            
            stats['chrom_counts'][chrom] += 1
            
            if filter_field != 'PASS':
                stats['filter_reasons']['non_PASS'] += 1
                stats['filtered_variants'] += 1
                continue
            
            format_keys = format_field.split(':')
            
            pass_samples = []
            for idx in tumor_indices:
                sample_values = sample_fields[idx].split(':')
                format_dict = dict(zip(format_keys, sample_values))
                
                if 'AD' not in format_dict:
                    stats['filter_reasons']['missing_AD'] += 1
                    continue
                
                ad_values = list(map(int, format_dict['AD'].split(',')))
                
                if len(ad_values) != len(alt.split(',')) + 1:
                    stats['filter_reasons']['unexpected_AD_count'] += 1
                    continue
                
                if any(ad >= min_alt_ad for ad in ad_values[1:]):
                    pass_samples.append(idx)
                    stats['sample_pass_counts'][samples[idx]] += 1
                    stats['chrom_sample_counts'][chrom][samples[idx]] += 1
            
            write_line = False
            if filter_mode == 'any' and pass_samples:
                write_line = True
            elif filter_mode == 'all' and len(pass_samples) == len(tumor_indices):
                write_line = True
            elif filter_mode == 'majority' and len(pass_samples) > len(tumor_indices) / 2:
                write_line = True
            
            if write_line:
                stats['passed_variants'] += 1
                pass_samples_info = ','.join([samples[i] for i in pass_samples])
                info_field = fields[7]
                info_field += f';PASS_SAMPLES={pass_samples_info}'
                fields[7] = info_field
                outfile.write('\t'.join(fields) + '\n')
            else:
                stats['filtered_variants'] += 1
                stats['filter_reasons']['AD_threshold'] += 1
    
    return stats

def write_log(log_file, stats, input_file, output_file, min_alt_ad, filter_mode):
    with open(log_file, 'w') as log:
        log.write(f"VCF Filtering Summary\n")
        log.write(f"=====================\n\n")
        log.write(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write(f"Input file: {os.path.basename(input_file)}\n")
        log.write(f"Output file: {os.path.basename(output_file)}\n")
        log.write(f"Minimum alternate allele depth: {min_alt_ad}\n")
        log.write(f"Filter mode: {filter_mode}\n")
        log.write(f"Filter mode explanation: {get_filter_mode_description(filter_mode)}\n")
        log.write(f"Processed samples: {', '.join(stats['processed_samples'])}\n\n")
        
        log.write(f"Variant Statistics:\n")
        log.write(f"  Total variants: {stats['total_variants']}\n")
        log.write(f"  Passed variants: {stats['passed_variants']}\n")
        log.write(f"  Filtered variants: {stats['filtered_variants']}\n")
        log.write(f"  Pass rate: {stats['passed_variants']/stats['total_variants']*100:.2f}%\n\n")
        
        log.write(f"Sample Pass Counts:\n")
        for sample in stats['processed_samples']:
            count = stats['sample_pass_counts'][sample]
            log.write(f"  {sample}: {count}\n")
        log.write("\n")
        
        log.write(f"Chromosome Counts:\n")
        log.write(f"{'Chromosome':<10}")
        for sample in stats['processed_samples']:
            log.write(f"{sample:<15}")
        log.write("Total\n")
        
        for chrom, total_count in stats['chrom_counts'].most_common():
            log.write(f"{chrom:<10}")
            for sample in stats['processed_samples']:
                count = stats['chrom_sample_counts'][chrom][sample]
                log.write(f"{count:<15}")
            log.write(f"{total_count}\n")
        log.write("\n")
        
        log.write(f"Filter Reasons:\n")
        for reason, count in stats['filter_reasons'].most_common():
            log.write(f"  {reason}: {count}\n")

def create_mutation_matrix(input_vcf, gene_list_file, output_tsv):
    # Read the list of genes
    with open(gene_list_file, 'r') as f:
        genes = [line.strip() for line in f]

    # Initialise the mutation matrix
    matrix = defaultdict(lambda: {gene: {} for gene in genes})
    
    # Process the VCF file
    with (gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, 'r')) as vcf:
        csq_format = None
        samples = []
        tumor_samples = []
        for line in vcf:
            if line.startswith('##tumor_sample='):
                tumor_samples.append(line.strip().split('=')[1])
            elif line.startswith('##INFO=<ID=CSQ'):
                csq_format = line.split('Format: ')[1].strip().strip('"').split('|')
            elif line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                break
        
        if not csq_format:
            raise ValueError("CSQ format not found in VCF header")
        
        if not tumor_samples:
            raise ValueError("No tumor samples found in VCF header")
        
        tumor_indices = [samples.index(sample) for sample in tumor_samples if sample in samples]
        if not tumor_indices:
            raise ValueError(f"No tumor samples found in VCF samples. Tumor samples in header: {', '.join(tumor_samples)}, VCF samples: {', '.join(samples)}")

        for line in vcf:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            info = dict(item.split('=') for item in fields[7].split(';') if '=' in item)
            
            if 'CSQ' in info:
                csq_entries = info['CSQ'].split(',')
                csq_data = dict(zip(csq_format, csq_entries[0].split('|')))  # Use the first CSQ annotation
                
                if csq_data['BIOTYPE'] == 'protein_coding' and csq_data['SYMBOL'] in genes:
                    gene = csq_data['SYMBOL']
                    consequence = csq_data['Consequence']
                    hgvsp = csq_data['HGVSp']
                    hgvsc = csq_data['HGVSc']
                    
                    important_consequences = {'missense_variant', 'stop_gained', 'frameshift_variant', 
                                              'splice_acceptor_variant', 'splice_donor_variant'}
                    if any(cons in consequence for cons in important_consequences):
                        if hgvsp:
                            change = hgvsp.split(':')[-1]
                        elif hgvsc:
                            change = hgvsc.split(':')[-1]
                        else:
                            change = 'Unknown'
                        
                        # Update the matrix for each tumor sample that has this variant
                        format_fields = fields[8].split(':')
                        gt_index = format_fields.index('GT')
                        af_index = format_fields.index('AF') if 'AF' in format_fields else None
                        
                        for i in tumor_indices:
                            sample_data = fields[9+i].split(':')
                            genotype = sample_data[gt_index]
                            if genotype not in ('0/0', '0|0', './.'):  # Check if the sample has the variant
                                af = sample_data[af_index] if af_index is not None else 'NA'
                                matrix[samples[i]][gene][change] = af

    # Write the matrix to a TSV file
    with open(output_tsv, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerow(['Sample'] + genes)
        
        for sample in tumor_samples:
            if sample in samples:  # Only include tumour samples that are actually in the VCF
                row = [sample]
                for gene in genes:
                    mutations = matrix[sample][gene]
                    if mutations:
                        mutation_str = ';'.join(f"{mut}({af})" for mut, af in mutations.items())
                        row.append(mutation_str)
                    else:
                        row.append('NA')
                writer.writerow(row)

def create_matrix_only(input_vcf, gene_list_file, output_tsv):
    print(f"Creating mutation matrix from {input_vcf}")
    create_mutation_matrix(input_vcf, gene_list_file, output_tsv)
    print(f"Mutation matrix saved to {output_tsv}")

def main():
    parser = argparse.ArgumentParser(description='Filter Mutect2 VCF files and/or create mutation matrix')
    parser.add_argument('input_file', help='Input VCF file (can be gzipped)')
    parser.add_argument('--output_file', help='Output VCF file for filtering (will be gzipped if ends with .gz)')
    parser.add_argument('--min_alt_ad', type=int, default=10, help='Minimum alternate allele depth (default: 10)')
    parser.add_argument('--all_samples', action='store_true', help='Process all samples instead of just the tumor samples')
    parser.add_argument('--sample', help='Process a specific sample by name')
    parser.add_argument('--filter_mode', choices=['any', 'all', 'majority'], default='any',
                        help='How to filter when multiple samples are present: '
                             'any (default, keep if any sample passes), '
                             'all (keep if all samples pass), '
                             'majority (keep if more than half of samples pass)')
    parser.add_argument('--log_file', help='Log file name (default: <output_file_basename>_filter_log.txt)')
    parser.add_argument('--gene_list', help='File containing list of genes to include in the matrix')
    parser.add_argument('--matrix_output', help='Output TSV file for the mutation matrix')
    parser.add_argument('--matrix_only', action='store_true', help='Create only the mutation matrix without filtering')
    
    args = parser.parse_args()
    
    if args.matrix_only:
        if not (args.gene_list and args.matrix_output):
            raise ValueError("Both --gene_list and --matrix_output must be provided when using --matrix_only")
        create_matrix_only(args.input_file, args.gene_list, args.matrix_output)
    else:
        if not args.output_file:
            raise ValueError("--output_file must be provided when not using --matrix_only")
        
        if args.all_samples and args.sample:
            raise ValueError("Cannot use both --all_samples and --sample options simultaneously")
        
        if not args.log_file:
            output_basename = os.path.splitext(os.path.basename(args.output_file))[0]
            args.log_file = f"{output_basename}_filter_log.txt"
        
        stats = filter_vcf(args.input_file, args.output_file, args.min_alt_ad, 
                           use_tumor_samples=(not args.all_samples and not args.sample),
                           specific_sample=args.sample,
                           filter_mode=args.filter_mode)
        
        write_log(args.log_file, stats, args.input_file, args.output_file, args.min_alt_ad, args.filter_mode)
        
        print(f"Filtered VCF saved to {args.output_file}")
        print(f"Filtering log saved to {args.log_file}")
        
        if args.gene_list and args.matrix_output:
            create_mutation_matrix(args.output_file, args.gene_list, args.matrix_output)
            print(f"Mutation matrix saved to {args.matrix_output}")
        elif args.gene_list or args.matrix_output:
            print("Warning: Both --gene_list and --matrix_output must be provided to create a mutation matrix.")

if __name__ == '__main__':
    main()