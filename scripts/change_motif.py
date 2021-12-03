# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import logging
import json, subprocess, re, os
os.chdir('/home/a7420174/tools/EHdn/scripts/')
from multiprocessing import Pool, cpu_count
from tempfile import NamedTemporaryFile
from core import regiontools
from functools import partial
from copy import copy
from multiprocessing import Pool, cpu_count

regions_path = '/home/a7420174/WGS/ASD_STR/outputs/ASD_EHdn.TR_list.id.tsv'
TRF_regions_path = '/home/a7420174/WGS/ASD_STR/resources/TRF_hg38_motif_N5.bed'
output_path = '/home/a7420174/WGS/ASD_STR/outputs/variant_catalog_hg38_custom.json'
info_path = '/home/a7420174/WGS/ASD_STR/outputs/EHdn_variant_info.json'


# %%
def load_ids(fname):
    ids = []
    with open(fname, "r") as bed_file:
        for line in bed_file:
            if not 'chr' in line:
                continue
            id = line.split()[-1]
            ids.append(id)
    return ids

def load_target_regions(fname):
    regions = []
    motifs = []
    with open(fname, "r") as bed_file:
        for line in bed_file:
            if not 'chr' in line:
                continue
            chrom, start, end, motif, *_ = line.split()
            start, end = int(start), int(end)
            region = regiontools.Region(chrom, start, end)
            regions.append(region)
            motifs.append(motif)
    return (regions, motifs)


# %%
trf_regions, trf_motifs = load_target_regions(TRF_regions_path)


# %%
trf_regions = load_ids(TRF_regions_path)

def MinimialUnitUnderShift(unit):
    minimal_unit = unit
    double_unit = unit + unit
    for index in range(0, len(unit)):
        current_unit = double_unit[index:index+len(unit)]
        if current_unit < minimal_unit:
            minimal_unit = current_unit
    return minimal_unit

def reverse_complement(unit):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    unit_rc = "".join(complement.get(base, base) for base in reversed(unit))
    return unit_rc

def ComputeCanonicalRepeatUnit(unit):
    minimal_unit = MinimialUnitUnderShift(unit)
    unit_rc = reverse_complement(unit)
    minimal_unit_rc = MinimialUnitUnderShift(unit_rc)

    if minimal_unit_rc < minimal_unit:
        return minimal_unit_rc
    return minimal_unit


# %%
with open(output_path, "r") as profile_file:
    catalog = json.load(profile_file)


# %%
def change_motif(Dict):
    region = Dict['ReferenceRegion']
    motif = Dict['LocusStructure'].strip('()*')
    for R, M in zip(trf_regions, trf_motifs):
        if (region == R) & (motif == MinimialUnitUnderShift(M)):
            Dict['LocusStructure'] = '(' + M + ')*'
    return Dict


# %%
with Pool(processes=40) as pool:
    new_catalog = pool.map(change_motif, catalog)


# %%
new_output_path = '/home/a7420174/WGS/ASD_STR/outputs/variant_catalog_hg38_custom2.json'


# %%
with open(new_output_path, 'w') as new_file:
    json.dump(new_catalog, new_file, indent='\t')


