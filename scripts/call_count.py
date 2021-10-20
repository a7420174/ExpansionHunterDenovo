#!/usr/bin/env python3

# Input: JSON file from ExpansionHunter Denovo
# Output: The number of raw calls (anchored in-repeat reasds) detected

import argparse

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("json_filename", type=str)
parser.add_argument("--ignore-alt-contigs", default=False, action="store_true")
parser.add_argument("--min-count", default=0, type=int)
args = parser.parse_args()
#####################################

def get_chromosomes():
    return([str(c) for c in list(range(1, 23)) + ["X", "Y"]])

chromosomes = get_chromosomes()

with open(args.json_filename) as f:
    read_depth = float(f.readlines()[-3].split(':')[1].strip().strip(','))

f = open(args.json_filename)
sample_count = 0


for line in f:
    if "RegionsWithIrrAnchors" in line:
        counting = True
    elif "RegionsWithIrrs" in line:
        counting = False
    elif "-" in line and counting:
        chrom = line.split(":")[0].replace('"', "").replace(" ", "")
        raw_count = float(line.split(':')[2].strip().strip(','))
        if not args.ignore_alt_contigs or chrom.replace("chr", "") in chromosomes:
            if (raw_count*40/read_depth) >= args.min_count:
                sample_count += 1
f.close()

print(sample_count)