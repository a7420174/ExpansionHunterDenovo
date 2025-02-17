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
from glob import glob
from tempfile import NamedTemporaryFile
from multiprocessing import Pool, cpu_count

number_threads = cpu_count() - 1
liftOver_path = '/home/a7420174/tools/UCSC/liftOver'
profile_list = glob('/home/a7420174/WGS/ASD_STR/inputs/ADNI_EHdn/*.str_profile.json')

def liftover_profile(profile):
    for _, rec in profile.items():
        if type(rec) != dict:
            continue
        if "RegionsWithIrrAnchors" not in rec:
            continue

        counts = []
        with NamedTemporaryFile('w', delete=False) as tmp:
            for region, count in rec["RegionsWithIrrAnchors"].items():
                tmp.write('chr'+region+'\n')
                counts.append(count)
            tmp_path = tmp.name

        new_region = []
        with NamedTemporaryFile('r') as tmp:
            tmp2_path = tmp.name
            cmd = ' '.join([
                    liftOver_path,
                    tmp_path,
                    '~/tools/UCSC/hg19ToHg38.over.chain.gz',
                    tmp2_path,
                    'unMapped',
                    '-positions'
            ])
            subprocess.run(cmd, shell=True)
            
            for line in tmp:
                new_region.append(line.strip())

        del rec["RegionsWithIrrAnchors"]
        rec["RegionsWithIrrAnchors"] = dict(zip(new_region, counts))
        os.remove(tmp_path)

    return profile


def run(profile_path):
    new_profile_path = re.sub('str_profile\.json', 'str_profile.hg38.json', profile_path)
    with open(profile_path, "r") as profile_file:
        profile = json.load(profile_file)
    new_profile = liftover_profile(profile)
    with open(new_profile_path, 'w') as new_file:
        json.dump(new_profile, new_file, indent='\t')

    logging.info("Done")

if __name__ == "__main__":
    pool = Pool(number_threads)
    pool.map(run, profile_list)
    pool.close()
    pool.join()
