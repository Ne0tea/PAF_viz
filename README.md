<<<<<<< HEAD
# PAF_viz

PAF_viz 是一个用于 **PAF 比对结果 dotplot 可视化** 的 Python 工具，支持区域限制、GFF3 基因结构注释、BED 三角标记与矩形高亮背景，适合基因组区段比对结果的快速展示。

## Features

- 支持 PAF 全局/局部区域 dotplot 绘制
- 支持 target/query 两侧 GFF3 基因结构轨道（gene/CDS/exon）
- 支持 BED 三角标记（triangle markers）与矩形背景高亮（rectangles）
- 支持按比对长度过滤，减少噪声
- 输出高质量 PDF 图（默认）
=======
<!--
 * @Descripttion: 
 * @Author: Ne0tea
 * @version: 
 * @Date: 2026-04-21 18:00:15
 * @LastEditors: Ne0tea
 * @LastEditTime: 2026-04-21 18:00:15
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
>>>>>>> d5e4d8b (feat: initial PAF_viz script with README and repo-scoped gitignore)

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
<<<<<<< HEAD
  -o all_query_target_dotplot
```

输出文件：

- `all_query_target_dotplot.pdf`

## Main Parameters

- `-i, --paf`: 输入 PAF 文件（必选）
- `-tr, --target-region`: target 区域，例如 `chr15:1000-5000`
- `-qr, --query-region`: query 区域，例如 `chr15:1000-5000`
- `-tg, -qg`: target/query GFF3 文件
- `-btx, -bty`: X/Y 轴 BED 三角标记文件
- `-brx, -bry`: X/Y 轴 BED 矩形背景文件
- `-flen`: 最小比对长度过滤阈值
- `-o`: 输出前缀
=======
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
>>>>>>> d5e4d8b (feat: initial PAF_viz script with README and repo-scoped gitignore)

## License

MIT
