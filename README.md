# PAF_viz

PAF_viz 是一个用于 **PAF 比对结果 dotplot 可视化** 的 Python 工具，支持区域限制、GFF3 基因结构注释、BED 三角标记与矩形高亮背景，适合基因组区段比对结果的快速展示。

## Features

- 支持 PAF 全局/局部区域 dotplot 绘制
- 支持 target/query 两侧 GFF3 基因结构轨道（gene/CDS/exon）
- 支持 BED 三角标记（triangle markers）与矩形背景高亮（rectangles）
- 支持按比对长度过滤，减少噪声
- 输出高质量 PDF 图（默认）

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

## License

MIT
