#!/usr/bin/env python3

import sys, argparse, logging

logging.basicConfig(level=logging.INFO)

# Initilize argument parser.
description = '''Sort allels in multi-allelic sites, such that the genotype is one of 0/0, 0/1 or 1/1.

Why? Because HapCUT2 exits with an error if it encounters a genotype such as 0/2. That is, genotypes containing
allele indicators larger than 1.

This script expects a single-sample VCF.

Say we have a variant such as the following:

chr1   1005        .       T       G,C     .       PASS     .  GT:AD:DP:PL       1/2:1,2,3:6:1,2,3,4,5,6

There are three alleles and the genotype is 1/2. The AD field contains one value per allele, and the PL field
contains one value per distinct genotype. The values in the PL field are ordered according to the genotypes in
the following order: 0/0, 0/1, 0/2, 1/1, 1/2, 2/2.

We want the genotype to be 0/1. In order to do this, we need the REF allele to be G and the ALT alleles to be C,T.
The result is the variant below.

chr1   1005        .       G       C,T     .       PASS     .  GT:AD:DP:PL       0/1:2,3,1:6:4,5,2,6,3,1

Note that we have changed the reference allele and the alternate allele order, the genotype, and the order of the
elements in the AD and PL fields.

HapCUT2 will not raise an error when it encounters the line above.
'''
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)

# Add the arguments.
parser.add_argument('--vcf', type=str, required=True,
        help='Path to single-sample VCF file. Not compressed.')

# Parse script input arguments.
args = parser.parse_args()

vcf = args.vcf


with open(vcf) as fid:
    for line in fid:
        if line[0] == '#':
            # Print the header lines just as they are.
            print(line, end='')
        else:
            # Remove trailing whitespace from line.
            line = line.strip()

            # Split the line into columns.
            columns = line.split('\t')

            # Get the reference allele.
            REF = column[3]

            # Get the list of alternate alleles.
            ALTS = column[4]
            # Split into a list.
            # Note that if there is only one alt allele, ALTS will be a list with a single element.
            ALTS = ALTS.split(',')

            if len(ALTS) < 2:
                # The site is diallelic, and does not need altering.
                # Print the unaltered line.
                print(line)

            # Make a list of alleles, where the first element is the reference allele, and the remaining
            # are the alternate alleles.
            ALLELES = [REF]
            ALLELES.extend(ALTS)

            # Get the FORMAT of the genotype field.
            format_column = columns[8]

            # Split the FORMAT column into individual fields which are separated by colons.
            format_fields = format_column.split(':')

            # The sample (or genotype) column.
            genotype_fields = columns[9]
            # Split into a list.
            genotype_fields = genotype_fields.split(':')

            # Get the sample genotype.
            GT_idx = format_fields.index('GT')
            GT = genotype_fields[GT_idx]

            # FIXME: I will definitely need to be smarter about parsing the genotype.
            # I can get inspiration from some of the other scripts.
            a1 = int(GT[0])
            a2 = int(GT[2])

            if min(a1, a2) < 2:
                # The genotype does not need altering, because it is either 0/0, 0/1 or 1/1.
                # Print the unaltered line.
                print(line)


            # Join the fields with colons again.
            new_format_column = ':'.join(format_fields)

            #new_row = ?????

            # Make the new line by joining all the columns by tabs.
            new_line = '\t'.join(new_row)

            print(new_line)

logging.warning('???')
