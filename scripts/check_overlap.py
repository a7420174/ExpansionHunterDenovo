# %%
import logging
import json, subprocess, re, os
from glob import glob
os.chdir('/home/a7420174/tools/EHdn/scripts/')
from multiprocessing import Pool, cpu_count
from tempfile import NamedTemporaryFile
from core import regiontools
from functools import partial
from copy import copy
from multiprocessing import Pool, cpu_count
import pandas as pd
import numpy as np

catalog_path = '/home/a7420174/WGS/ASD_STR/inputs/EH_1kg_row_table.tsv'
catalog2_path = '/home/a7420174/WGS/ASD_STR/inputs/EH_1kg_custom_row_table.tsv'
json_list = glob('/home/a7420174/WGS/ASD_STR/inputs/1kg_EH_custom/*.json')
vcf_list = glob('/home/a7420174/WGS/ASD_STR/inputs/1kg_EH_custom/*.vcf')
number_threads = cpu_count()

# %%
def load_target_regions(fname):
    regions = []
    motifs = []
    ids = []
    with open(fname, "r") as bed_file:
        for line in bed_file:
            if not 'chr' in line:
                continue
            locus_id, var_id, region, motif = line.strip().split('\t')
            chrom, interval = region.split(':')
            start, end = interval.split('-')
            id = '_'.join([chrom, start, end, motif])
            start, end = int(start), int(end)
            region = regiontools.Region(chrom, start, end)
            regions.append(region)
            motifs.append(motif)
            ids.append(id)
    return (regions, motifs, ids)

# %%
regions, motifs, ids = load_target_regions(catalog_path)
regions2, motifs2, ids2 = load_target_regions(catalog2_path)

# %%
def overlap(region, motif, regions, motifs, ids):
    overlap_id = []
    for R, ID, M in zip(regions, ids, motifs):
        if (regiontools.compute_distance(region, R) == 0):# & (motif == M):
            overlap_id.append(ID)
    overlap_id = ':'.join(overlap_id)
    return overlap_id

#%%
overlap_partial = partial(overlap, regions = regions2, ids = ids2, motifs = motifs2)
with Pool(processes=number_threads) as pool:
    overlap_var_ids = pool.starmap(overlap_partial, zip(regions, motifs))
#%%
compar = pd.DataFrame({'illumina': ids, 'custom': overlap_var_ids})
compar = compar[compar.custom != '']
compar['custom'] = compar.custom.apply(lambda x: x.split(':'))
#%%
flatdata = pd.DataFrame([(index, value) for (index, values)
                        in compar['custom'].iteritems() for value in values],
                            columns = ['index', 'custom']).set_index( 'index' )

compar = compar.drop('custom', axis = 1 ).join( flatdata )
#%%
def MinimialUnitUnderShift(unit):
    minimal_unit = unit
    double_unit = unit + unit
    for index in range(0, len(unit)):
        current_unit = double_unit[index:index+len(unit)]
        if current_unit < minimal_unit:
            minimal_unit = current_unit
    return minimal_unit

compar['illumina_motif'] = compar.illumina.str.replace('.*_([ATGCNR]+)$', '\\1').apply(MinimialUnitUnderShift)
compar['custom_motif'] = compar.custom.str.replace('.*_([ATGCNR]+)$', '\\1').apply(MinimialUnitUnderShift)
overlap_locus = compar[compar.illumina_motif == compar.custom_motif].custom.tolist()
exceptions = ['chr4_39348424_39348483_AAAAG', 'chr4_39348424_39348483_AAGGG', 'chr4_41745975_41746037_GCC']
overlap_locus.extend(exceptions)

#%%
idx = []
for i, id in enumerate(ids2):
    if id in overlap_locus:
        idx.append(i)
# %%
catalog = pd.read_csv(catalog_path, sep='\t')
catalog2 = pd.read_csv(catalog2_path, sep='\t')
# %%
new_catalog2 = catalog2.iloc[~catalog2.index.isin(idx)]
concat_catalog= pd.concat([catalog, new_catalog2]).reset_index(drop=True)
concat_catalog[['contig', 'start', 'end']]  = concat_catalog.region.str.split('[:-]', expand=True)
concat_catalog = concat_catalog.astype({'contig': 'category', 'start': 'int', 'end': 'int'})
concat_catalog['contig'] = concat_catalog.contig.cat.reorder_categories(['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY'], ordered=True)
concat_catalog = concat_catalog.sort_values(['contig', 'start', 'end']).drop(columns=['contig', 'start', 'end'])
var_idx = concat_catalog.index.tolist()

# %%
concat_catalog.to_csv('/home/a7420174/WGS/ASD_STR/inputs/EH_KOR_ASD_row_table_final.tsv', sep='\t', index=False)
concat_catalog.to_csv('/home/a7420174/WGS/ASD_STR/inputs/EH_1kg_row_table_final.tsv', sep='\t', index=False)
#%%
cohort = '1kg'

GT1 = np.genfromtxt('/home/a7420174/WGS/ASD_STR/inputs/EH_'+cohort+'_GT_matrix.tsv.gz',
              dtype='str', delimiter='\t', missing_values='NA')

GT2 = np.genfromtxt('/home/a7420174/WGS/ASD_STR/inputs/EH_'+cohort+'_custom_GT_matrix.tsv.gz',
              dtype='str', delimiter='\t', missing_values='NA')
#%%
GT = np.concatenate((GT1, GT2), axis=0)
#%%
new_GT = GT[var_idx,]
#%%
np.savetxt('/home/a7420174/WGS/ASD_STR/inputs/EH_'+cohort+'_GT_matrix_final.tsv.gz', new_GT, delimiter='\t', fmt='%s')

# %%
# duplicated = []
# for id in set(overlap_var_ids):
#     if re.search(':', id):
#         if len(id.split(':')) > 1:
#             overlap_id = id.split(':')[-1]
#             duplicated.append(overlap_id)
#             if len(id.split(':')) > 2:
#                 overlap_id = id.split(':')[-2]
#                 duplicated.append(overlap_id)
#                 if len(id.split(':')) > 3:
#                     overlap_id = id.split(':')[-3]
#                     duplicated.append(overlap_id)
#%%
# def mod_json(profile_path):
#     new_profile_path = re.sub('\.json', '.mod.json', profile_path)
#     with open(profile_path, "r") as profile_file:
#         profile = json.load(profile_file)
#     for id in duplicated:
#         del profile['LocusResults'][id]
#     with open(new_profile_path, 'w') as new_file:
#         json.dump(profile, new_file, indent='\t')
# #%%
# with Pool(processes=number_threads) as pool:
#     pool.map(mod_json, json_list)
# # %%
# def mod_vcf(vcf_path):
#     new_vcf_path = re.sub('\.vcf', '.mod.vcf', vcf_path)
#     new_lines = []
#     with open(vcf_path, "r") as vcf:
#         for line in vcf:
#             if re.match('#', line):
#                 new_lines.append(line)
#             else:
#                 var_id = re.findall('VARID=(chr[0-9XY]+_[0-9]+_[0-9]+_[ATGC]+)', line)[0]
#                 if var_id not in duplicated:
#                         new_lines.append(line)
#     with open(new_vcf_path, "w") as new_vcf:
#         new_vcf.writelines(new_lines)
# #%%
# with Pool(processes=number_threads) as pool:
#     pool.map(mod_vcf, vcf_list)
