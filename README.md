<!--
 * @Descripttion: 
 * @Author: Ne0tea
 * @version: 
 * @Date: 2026-04-21 18:00:15
 * @LastEditors: Ne0tea
 * @LastEditTime: 2026-04-21 18:19:08
-->
# PAF_viz

PAF_viz is a Python tool for visualizing **PAF alignment results as dot plots**.  
It supports region-specific plotting, GFF3-based gene structure tracks, BED triangle markers, and BED rectangle highlights, making it suitable for quick inspection of genome segment alignments.

## Features

- Plot global or region-restricted dot plots from PAF files
- Display gene structure annotation tracks (gene/CDS/exon) from GFF3 on both target and query axes
- Add BED-based triangle markers on axes
- Add BED-based rectangle highlight backgrounds
- Filter short alignments by minimum alignment length to reduce noise
- Export high-quality figures (PDF by default)

## Requirements

- Python 3.8+
- matplotlib
- pandas

## Installation

```bash
pip install matplotlib pandas
```

## Quick Start

```bash
python PAF_viz.py \
  -i your_alignment.paf \
  -tg target.gff3 \
  -qg query.gff3 \
  -btx query_marks.bed \
  -bty target_marks.bed \
  -brx query_rects.bed \
  -bry target_rects.bed \
  -tr chr15:17186630-31196487 \
  -flen 2000 \
  -o PREFIX
```

Output file:

- `PREFIX.pdf`

## Main Arguments

| Argument | Description |
| --- | --- |
| `-i, --paf` | Input PAF file (**required**) |
| `-tr, --target-region` | Target region, e.g. `chr15:1000-5000` |
| `-qr, --query-region` | Query region, e.g. `chr15:1000-5000` |
| `-tg` | Target-side GFF3 annotation file |
| `-qg` | Query-side GFF3 annotation file |
| `-btx` | BED file for X-axis triangle markers |
| `-bty` | BED file for Y-axis triangle markers |
| `-brx` | BED file for X-axis rectangle highlights |
| `-bry` | BED file for Y-axis rectangle highlights |
| `-flen` | Minimum alignment length threshold |
| `-o` | Output prefix |

## Input File Types

- **PAF**: pairwise alignment records (primary plotting source)
- **GFF3** *(optional)*: gene structure annotations
- **BED** *(optional)*: axis markers and highlight intervals
