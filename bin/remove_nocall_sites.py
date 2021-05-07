#!/usr/bin/env python3

import sys, argparse, logging

logging.basicConfig(level=logging.INFO)

# Initilize argument parser.
description = '''Remove variants with no-call genotypes in single sample VCFs.

Sites where the genotype is "./." will be removed.
'''
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)

# Add the arguments.
parser.add_argument('--vcf', type=str, required=True,
        help='Path to VCF file. Not compressed.')

# Parse script input arguments.
args = parser.parse_args()

vcf = args.vcf

with open(vcf) as fid:
    for line in fid:
        if line[:2] == '##':
            print(line, end='')
            
        else:
            # Remove trailing whitespace from line.
            line = line.strip()

            # Split the line into columns.
            columns = line.split('\t')

            # The sample (or genotype) columns.
            sample_columns = columns[9:]
            
            # If the current line is the table header, we will check that the VCF is single sample.
            if line[:6] == '#CHROM':
                assert len(sample_columns) == 1, 'This script only accepts single sample VCFs. This VCF contains {n} samples.'.format(n=len(sample_columns))
            
            # Split the values.
            sample_fields = sample_columns[0].split(':')
            # Get the genotype value.
            GT = sample_fields[0]
            
            # If the variant is not called in the sample, it is coded as "./.".
            if GT == './.':
                continue
            
            print(line)
