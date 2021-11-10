#
# ExpansionHunter Denovo
# Copyright 2016-2019 Illumina, Inc.
# All rights reserved.
#
# Author: Egor Dolzhenko <edolzhenko@illumina.com>,
#         Michael Eberle <meberle@illumina.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#

import logging
import json, subprocess, re, os
from multiprocessing import Pool, cpu_count
from tempfile import NamedTemporaryFile
from core import regiontools
from functools import partial
from copy import copy


regions_path = '/home/a7420174/WGS/ASD_STR/outputs/ASD_EHdn.TR_list.id.tsv'
TRF_regions_path = '/home/a7420174/WGS/ASD_STR/resources/TRF_hg38_motif.bed'
output_path = '/home/a7420174/WGS/ASD_STR/outputs/variant_catalog_hg38_custom.json'
info_path = '/home/a7420174/WGS/ASD_STR/outputs/EHdn_variant_info.json'

number_threads = cpu_count()

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

def load_ids(fname):
    ids = []
    with open(fname, "r") as bed_file:
        for line in bed_file:
            if not 'chr' in line:
                continue
            id = line.split()[-1]
            ids.append(id)
    return ids

def MinimialUnitUnderShift(unit):
    minimal_unit = unit
    double_unit = unit + unit
    for index in range(0, len(unit)):
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
        return minimal_unit_rc, 1
    return minimal_unit, 0

def match_regions(region, motif, id, regions, motifs):
    matched_regions = []
    for R, M in zip(regions, motifs):
        canonical_M, is_rc = ComputeCanonicalRepeatUnit(M)
        Distance = regiontools.compute_distance(region, R)
        Cond = (Distance <= 1000) & (motif in canonical_M)
        if Cond:
            matched_region = R.chrom + ':' + str(R.start) + '-' + str(R.end)
            matched_regions.append([matched_region, canonical_M, is_rc, Distance])
    if len(matched_regions) == 0:
        return id, 'NA', motif, 'NA'
    Dist = 1000
    for matched_region, canonical_M, is_rc, Distance in matched_regions:
        if motif == canonical_M:
            if is_rc == 1:
                return id, matched_region, reverse_complement(motif), 'Equal'
            else:
                return id, matched_region, motif, 'Equal'
        if Distance < Dist:
            Dist = Distance
    for matched_region, canonical_M, is_rc, Distance in matched_regions:
        if Distance == Dist:
            if is_rc == 1:
                return id, matched_region, reverse_complement(motif), 'Included'
            else:
                return id, matched_region, motif, 'Included'


def main():
    regions, motifs = load_target_regions(regions_path)
    TRF_regions, TRF_motifs = load_target_regions(TRF_regions_path)
    ids = load_ids(regions_path)
    match_regions_partial = partial(match_regions, regions = TRF_regions, motifs = TRF_motifs)
    with Pool(processes=number_threads) as pool:
        matched_regions = pool.starmap(match_regions_partial, zip(regions, motifs, ids))

    catalog = []
    for id, region, motif, _ in matched_regions:
        if region != 'NA':
            json_object = dict()
            json_object['LocusId'] = id
            json_object['LocusStructure'] = '(' + motif + ')*'
            json_object['ReferenceRegion'] = region
            json_object['VariantType'] = 'Repeat'
            catalog.append(json_object)


    info = []
    for id, region, motif, status in matched_regions:
        json_object = dict()
        json_object['LocusId'] = id
        json_object['Motif'] = motif
        json_object['ReferenceRegion'] = region
        json_object['Feature'] = status
        info.append(json_object)

    with open(info_path, 'w') as new_file:
        json.dump(info, new_file, indent='\t')

    logging.info("Done")
    return 0


if __name__ == "__main__":
    main()

