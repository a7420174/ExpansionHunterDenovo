# ExpansionHunter Denovo (Custom)

ExpansionHunter Denovo (EHdn) is a suite of tools for detecting novel expansions
of short tandem repeats (STRs). EHdn is intended for analysis of a collection of
BAM/CRAM files containing alignments of short (100-200bp) reads.

## 기존 EHdn v0.9와 다른 점

### 스크립트 추가 및 수정

1. list_TR.py
  
  Tandem repeat 목록을 생성하는 스크립트

- Column info: Tandem repeat의 locus, motif, motif의 길이와 A,T nucleotide의 비율, Anchored IIRs에 대한 모든 샘플의 count의 합

- Usage

```bash
list_TR.py locus \
--multisample-profile multisample_profile.json \
--output TR_list.tsv \            
--min-count 2
```

2. outlier.py

- `--min-count`, `--z-score`, `--no-ctrl` option 추가 (locus 분석만 가능하고 motif 분석은 안됨)

  - `--min-count`: The minimum of depth-normalized count in at least one sample (Default: 2)

  - `--z-score`: The cutoff of z score (Default: 1)

  - `--no-ctrl`: If at least one control has z score higher than cutoff, skip.
  options:  0(false) or 1(true) (Default: 0)

- TSV 파일 형식 변환

  - `counts` column에서 count에 대한 sample id가 나오도록 변경

  - `top_case_zscore` column 대신 `high_case_zscores` column에 cutoff를 넘는 모든 sample의 id와 z score 표시

- Usage

```bash
outlier.py locus --manifest manifest.txt \
--multisample-profile multisample_profile.json \
--output outlier_locus.tsv \
--min-count 2 \
--z-score 10 \
--no-ctrl 1
```

3. hg19_to_hg38.py

- `liftOver_path`: The path of UCSC liftOver
- `profile_list`: List of EHdn STR profiles

## License

ExpansionHunter Denovo is provided under the terms and conditions of the [Apache
License Version 2.0](LICENSE.txt). It relies on several third party packages
provided under other open-source licenses, please see [COPYRIGHT.txt](COPYRIGHT.txt)
for additional details.
