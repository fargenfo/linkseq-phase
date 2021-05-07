#!/usr/bin/env python3

import sys, argparse, logging

logging.basicConfig(level=logging.INFO)

# Initilize argument parser.
description = '''Remove a field from the FORMAT of a VCF file.

For example, if you remove "DP", the FORMAT line in the header with ID=DP will be removed, and
the FORMAT of all variants with this field will be updated.

Say the FORMAT of a variant is as follows: GT:DP:AD. Then the new FORMAT is GT:AD.
If a sample has genotype 0/1:40:18, then the new value will be 0/1:18.

Only one field may be supplied
'''
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)

# Add the arguments.
parser.add_argument('--vcf', type=str, required=True,
        help='Path to VCF file. Not compressed.')
parser.add_argument('--field', type=str, required=True,
        help='Name of field to remove from format. For example "DP".')

# Parse script input arguments.
args = parser.parse_args()

vcf = args.vcf
field = args.field

format_line_removed = False

with open(vcf) as fid:
    for line in fid:
        if line[0] == '#':
            # The header will contain a line like this:
            ##FORMAT=<ID=[FIELD], [...]>
            # We will remove a "FORMAT" line with "ID" matching the field we want to remove.
            if line[:8] == '##FORMAT' and 'ID=' + field in line:
                format_line_removed = True
                continue

            print(line, end='')
        else:
            # Remove trailing whitespace from line.
            line = line.strip()

            # Split the line into columns.
            columns = line.split('\t')

            # The 8th column is the FORMAT column.
            format_column = columns[8]

            # The sample (or genotype) columns.
            sample_columns = columns[9:]

            # Start building the new row.
            new_row = columns[:8]

            if field in format_column:
                # Split the FORMAT column into individual fields which are separated by colons.
                format_fields = format_column.split(':')

                # Get the index of the relevant field in the list.
                field_idx = format_fields.index(field)

                # Remove this field from the list of fields.
                format_fields.pop(field_idx)

                # Join the fields with colons again.
                new_format_column = ':'.join(format_fields)

                # Start making the new line.
                new_row.append(new_format_column)

                # Remove this field from each of the samples genotype fields.
                for sample_column in sample_columns:
                    # Split the column by colon into fields.
                    sample_fields = sample_column.split(':')
                    # Remove the field.
                    sample_fields.pop(field_idx)
                    # Join the fields back with colons.
                    new_sample_fields = ':'.join(sample_fields)
                    # Append the column to the other columns.
                    new_row.append(new_sample_fields)

                # Make the new line by joining all the columns by tabs.
                new_line = '\t'.join(new_row)
            else:
                # If the relevant field is not present in the line, just print the line.
                new_line = line

            print(new_line)

if not format_line_removed:
    logging.warning('A FORMAT line matching the {field} field was not found in the header of the VCF'.format(field=field))
