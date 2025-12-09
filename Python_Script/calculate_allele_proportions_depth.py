import argparse
import re
import sys

def calculate_proportions(vcf_file_path, output_tsv_path, summary_file_path, min_depth):
    """
    Parses a single-sample VCF file, filters sites based on depth and reference
    read count, calculates the proportion of alternative reads, and writes the
    results to output files.

    Args:
        vcf_file_path (str): The path to the input VCF file.
        output_tsv_path (str): The path for the output per-position TSV file.
        summary_file_path (str): The path for the output summary statistics file.
        min_depth (int): The minimum read depth (DP) to filter sites.
    """
    line_number = 0
    records_processed = 0
    sites_skipped_depth = 0
    sites_skipped_ref_zero = 0
    
    # Variables for genome-wide summary
    genome_total_reads = 0
    genome_total_alt_reads = 0

    print(f"Processing VCF file: {vcf_file_path}")
    print(f"Using minimum read depth filter: {min_depth}")
    print(f"Per-position output will be saved to: {output_tsv_path}")
    print(f"Summary output will be saved to: {summary_file_path}")

    try:
        with open(vcf_file_path, 'r') as vcf_in, open(output_tsv_path, 'w') as tsv_out:
            # Write the header for the per-position TSV file
            header = "CHROM\tPOS\tID\tREF\tALT\tTotalReads\tAltReads\tProportion\n"
            tsv_out.write(header)

            for line in vcf_in:
                line_number += 1
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 10:
                    print(f"Warning: Skipping malformed line {line_number}. Expected at least 10 columns.", file=sys.stderr)
                    continue

                # --- 1. Extract basic VCF columns ---
                chrom, pos, rs_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
                info_field = fields[7]
                format_field = fields[8]
                sample_field = fields[9]

                total_reads = 0
                alt_reads = 0
                ref_reads = 0

                # --- 2. Get Total Read Count and apply depth filter ---
                dp_match = re.search(r'DP=(\d+)', info_field)
                if dp_match:
                    total_reads = int(dp_match.group(1))
                    # FILTER 1: Skip sites with read depth < min_depth
                    if total_reads < min_depth:
                        sites_skipped_depth += 1
                        continue
                else:
                    print(f"Warning: DP tag not found on line {line_number}. Skipping site.", file=sys.stderr)
                    continue

                # --- 3. Get Allele Depths (AD) and apply zero-reference filter ---
                format_parts = format_field.split(':')
                try:
                    ad_index = format_parts.index('AD')
                    sample_parts = sample_field.split(':')
                    
                    if ad_index < len(sample_parts):
                        ad_values_str = sample_parts[ad_index]
                        if ad_values_str != '.':
                            ad_counts = [int(count) for count in ad_values_str.split(',')]
                            
                            # FILTER 2: Skip sites with zero reference reads
                            ref_reads = ad_counts[0]
                            if ref_reads == 0:
                                sites_skipped_ref_zero += 1
                                continue

                            if len(ad_counts) > 1:
                                alt_reads = sum(ad_counts[1:])
                    else:
                         print(f"Warning: Mismatch between FORMAT and sample fields on line {line_number}. Cannot find AD value.", file=sys.stderr)
                         continue

                except ValueError:
                    print(f"Warning: AD tag not found in FORMAT field on line {line_number}. Skipping site.", file=sys.stderr)
                    continue
                
                # --- 4. Calculate Proportion ---
                # Recalculate total_reads from AD to be precise, as DP can sometimes differ slightly.
                effective_total_reads = ref_reads + alt_reads
                
                # Double-check effective depth against the filter, as AD might sum to less than DP
                if effective_total_reads < min_depth:
                     sites_skipped_depth += 1
                     continue

                proportion = 0.0
                if effective_total_reads > 0:
                    proportion = alt_reads / effective_total_reads
                
                # --- 5. Write to per-position output file ---
                output_line = f"{chrom}\t{pos}\t{rs_id}\t{ref}\t{alt}\t{effective_total_reads}\t{alt_reads}\t{proportion:.4f}\n"
                tsv_out.write(output_line)
                records_processed += 1

                # --- 6. Update genome-wide totals ---
                genome_total_reads += effective_total_reads
                genome_total_alt_reads += alt_reads

    except FileNotFoundError:
        print(f"Error: The file '{vcf_file_path}' was not found.", file=sys.stderr)
        return
    except Exception as e:
        print(f"An unexpected error occurred while processing: {e}", file=sys.stderr)
        return

    # --- 7. Write the summary file ---
    try:
        with open(summary_file_path, 'w') as summary_out:
            summary_out.write("Total_Reads\tTotal_Alternative_Reads\tRatio_Alt_vs_Total\n")
            
            genome_ratio = 0.0
            if genome_total_reads > 0:
                genome_ratio = genome_total_alt_reads / genome_total_reads
            
            summary_out.write(f"{genome_total_reads}\t{genome_total_alt_reads}\t{genome_ratio:.6f}\n")
    except Exception as e:
        print(f"An error occurred while writing the summary file: {e}", file=sys.stderr)

    print("\n--- Processing Summary ---")
    print(f"Sites skipped due to read depth < {min_depth}: {sites_skipped_depth}")
    print(f"Sites skipped due to zero reference reads: {sites_skipped_ref_zero}")
    print(f"Total variant records written to file: {records_processed}")
    print("--------------------------")

def main():
    """
    Main function to handle command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Filter and calculate allele proportions from a single-sample VCF file."
    )
    # Positional arguments
    parser.add_argument("vcf_file", help="Path to the input VCF file.")
    parser.add_argument("output_tsv", help="Path for the output per-position TSV file.")
    parser.add_argument("summary_file", help="Path for the output summary statistics file.")
    
    # Optional argument for min depth
    parser.add_argument(
        "-d", "--min-depth",
        type=int,
        default=10,
        help="Minimum read depth (DP) required to process a site. (default: 10)"
    )
    
    args = parser.parse_args()

    calculate_proportions(args.vcf_file, args.output_tsv, args.summary_file, args.min_depth)

if __name__ == "__main__":
    main()
