'''
Descripttion: PAF alignment dotplot visualization with enhanced features
Author: Ne0tea
version: 3.0
Date: 2024-05-27 19:21:44
LastEditors: Ne0tea
LastEditTime: 2026-04-24 16:08:27
'''
import argparse
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import logging
import math
import re

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import transforms as mtransforms
from matplotlib.colors import is_color_like
from matplotlib.patches import Rectangle

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# Default colors
RECT_ALPHA = 0.3
TITLE_SIZE = 12
LABEL_SIZE = 7
GENE_SIZE = 7
TRIANGLE_MARKER_SIZE_PT = 6
TRIANGLE_TRACK_PAD_PT = 0.5
GENE_TRACK_PITCH = 1.5
GENE_TRACK_SLOT_PT = 7.0
GENE_TRACK_PAD_PT = 1.5

COLORS = {
    'triangle_default': '#780000',
    'rect_default': '#0066CC',
    'gene_cds': '#666666',
    'gene_exon': '#888888',
}
LAYOUT = {
    'delta_factor': 0.005,
    'label_pad_factor': 0.003,
    'track_label_pad': 0.18,
    'max_feat_span': 0.55,
    'track_margin_min': 0.04,
}

Region = Optional[Tuple[str, int, int]]


@dataclass
class PlotInputConfig:
    """CLI 输入与运行参数配置。"""
    paf_file: str
    target_region: Optional[object] = None
    query_region: Optional[object] = None
    target_gff: Optional[str] = None
    query_gff: Optional[str] = None
    bed_tri_x: Optional[str] = None
    bed_tri_y: Optional[str] = None
    bed_rect_x: Optional[str] = None
    bed_rect_y: Optional[str] = None
    filter_len: int = 2000
    output_prefix: str = 'all_query_target_dotplot'


@dataclass
class BedCollections:
    tri_x: List[Tuple[str, int, int, Optional[str]]] = field(default_factory=list)
    tri_y: List[Tuple[str, int, int, Optional[str]]] = field(default_factory=list)
    rect_x: List[Tuple[str, int, int, Optional[str]]] = field(default_factory=list)
    rect_y: List[Tuple[str, int, int, Optional[str]]] = field(default_factory=list)

    @classmethod
    def from_files(
        cls,
        bed_tri_x: Optional[str],
        bed_tri_y: Optional[str],
        bed_rect_x: Optional[str],
        bed_rect_y: Optional[str],
    ) -> 'BedCollections':
        return cls(
            tri_x=parse_bed_file(bed_tri_x),
            tri_y=parse_bed_file(bed_tri_y),
            rect_x=parse_bed_file(bed_rect_x),
            rect_y=parse_bed_file(bed_rect_y),
        )


@dataclass
class AxisLayout:
    names: List[str]
    lengths: Dict[str, int]
    offsets: Dict[str, int]
    total: int
    region: Region = None

    @classmethod
    def build(
        cls,
        dataframe,
        name_col: str,
        length_col: str,
        region: Optional[object] = None,
    ) -> 'AxisLayout':
        normalized_region = normalize_region(region)
        names, lengths, offsets, total = build_axis_layout(
            dataframe,
            name_col,
            length_col,
            region=normalized_region,
        )
        return cls(
            names=names,
            lengths=lengths,
            offsets=offsets,
            total=total,
            region=normalized_region,
        )


@dataclass
class GeneDrawContext:
    seq_offsets: Dict[str, int]
    axis_scale: float
    region: Region = None


@dataclass
class PlotRuntimeContext:
    dataframe: pd.DataFrame
    target_layout: AxisLayout
    query_layout: AxisLayout
    target_genes: List[dict]
    query_genes: List[dict]
    beds: BedCollections
    axis_scale: float
    unit_label: str
    output_prefix: str

def _safe_color(color_value, default_color):
    """Return a valid matplotlib color, otherwise fallback to default."""
    if color_value is None:
        return default_color
    c = str(color_value).strip()
    if not c:
        return default_color
    if is_color_like(c):
        return c
    return default_color

def _merge_gene_features(cds_list, exon_list, gene_list):
    """Construct gene structure from CDS/exon."""
    if gene_list:
        return gene_list
    
    merged = {}
    for feat in cds_list + exon_list:
        gid = feat['gene_id']
        if not gid:
            continue
        if gid not in merged:
            merged[gid] = {
                'chrom': feat['chrom'],
                'start': feat['start'],
                'end': feat['end'],
                'gene_id': gid,
            }
        else:
            merged[gid]['start'] = min(merged[gid]['start'], feat['start'])
            merged[gid]['end'] = max(merged[gid]['end'], feat['end'])
    return list(merged.values())

def _filter_gene_features(genes):
    """Split gene features."""
    cds_list = [g for g in genes if g['type'] == 'CDS' and g['gene_id']]
    exon_list = [g for g in genes if g['type'] == 'exon' and g['gene_id']]
    gene_list = [g for g in genes if g['type'] == 'gene']
    return cds_list, exon_list, gene_list


def parse_region(region_str):
    """Parse region string format: chrom:start-end"""
    match = re.match(r'([^:]+):(\d+)-(\d+)', region_str)
    if match:
        chrom, start, end = match.groups()
        return chrom, int(start), int(end)
    raise ValueError(f"Invalid region format: {region_str}. Expected format: chrom:start-end")


def normalize_region(region):
    """Normalize region to (chrom, start, end) with start < end.

    Accepts string: chrom:start-end or tuple/list: (chrom, start, end).
    """
    if region is None:
        return None

    if isinstance(region, str):
        chrom, start, end = parse_region(region)
    elif isinstance(region, (tuple, list)) and len(region) == 3:
        chrom, start, end = region
    else:
        raise ValueError(
            f'Invalid region value: {region}. Expected "chrom:start-end" or (chrom, start, end).'
        )

    start = int(start)
    end = int(end)
    if start == end:
        raise ValueError(f'Invalid region {chrom}:{start}-{end}: start must not equal end.')

    start, end = sorted((start, end))
    return str(chrom), start, end


def ensure_region_exists_in_paf(dataframe, region, axis='target'):
    """Validate region exists in PAF for axis (target/query), return normalized region."""
    normalized = normalize_region(region)
    if normalized is None:
        return None

    chrom, start, end = normalized
    name_col = f'{axis}_name'
    start_col = f'{axis}_start'
    end_col = f'{axis}_end'
    length_col = f'{axis}_length'

    chrom_df = dataframe[dataframe[name_col].astype(str) == chrom]
    if chrom_df.empty:
        raise ValueError(f'{axis}_region {chrom}:{start}-{end} does not exist in PAF ({name_col} not found).')

    seq_len = int(chrom_df[length_col].max())
    clipped_start = max(0, min(start, seq_len))
    clipped_end = max(0, min(end, seq_len))
    if clipped_end <= clipped_start:
        raise ValueError(
            f'{axis}_region {chrom}:{start}-{end} is outside sequence range [0, {seq_len}] in PAF.'
        )

    overlap_mask = (chrom_df[end_col] >= clipped_start) & (chrom_df[start_col] <= clipped_end)
    if not overlap_mask.any():
        raise ValueError(
            f'{axis}_region {chrom}:{start}-{end} has no overlapping alignments in PAF.'
        )

    return chrom, clipped_start, clipped_end


def parse_bed_file(bed_file):
    """Parse BED file and return list of (chrom, start, end, color) tuples"""
    bed_data = []
    if bed_file is None:
        return bed_data

    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            color = parts[3] if len(parts) >= 4 else None
            bed_data.append((chrom, start, end, color))

    return bed_data


def parse_gff3_attributes(attr_str):
    """Parse GFF3 attributes string into dict"""
    attrs = {}
    for item in attr_str.split(';'):
        item = item.strip()
        if not item or '=' not in item:
            continue
        key, value = item.split('=', 1)
        attrs[key] = value
    return attrs


def parse_gff3(gff_file):
    """Parse GFF3 file and resolve gene/CDS/exon relationship by ID/Parent"""
    if gff_file is None:
        return []

    features = []
    with open(gff_file, 'r') as f:
        for line in f:
            raw_line = line.rstrip('\n')
            if not raw_line.strip():
                continue
            if raw_line.lstrip().startswith('#'):
                continue

            # Rule for gene ID labels:
            # only when a GFF3 feature line is explicitly marked by leading '*'
            # (e.g. "*chr...\t...\tgene\t...") should that gene show ID label.
            line_show_gene_id = False
            line_for_parse = raw_line.strip()
            if line_for_parse.startswith('*'):
                line_show_gene_id = True
                line_for_parse = line_for_parse[1:].lstrip()

            parts = line_for_parse.split('\t')
            if len(parts) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attr_text = parts[:9]
            if feature_type.startswith('*'):
                line_show_gene_id = True
                feature_type = feature_type.lstrip('*')

            # Strict rule for feature color:
            # only read color from the last extra column (10th+), not from attributes.
            cur_color = parts[-1].strip() if len(parts) >= 10 else None
            attrs = parse_gff3_attributes(attr_text)
            feature_id = attrs.get('ID')
            parents = attrs.get('Parent', '')
            parent_ids = [p for p in parents.split(',') if p] if parents else []

            features.append({
                'chrom': chrom,
                'type': feature_type,
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'id': feature_id,
                'parents': parent_ids,
                'gene_id': None,
                'color': cur_color,
                'show_gene_id': bool(line_show_gene_id and feature_type == 'gene'),
            })

    if not features:
        return features

    gene_ids = {f['id'] for f in features if f['type'] == 'gene' and f['id']}
    transcript_types = {
        'mRNA', 'transcript', 'lnc_RNA', 'ncRNA', 'tRNA',
        'rRNA', 'miRNA', 'snRNA', 'snoRNA', 'pseudogenic_transcript'
    }

    transcript_to_gene = {}

    # First pass: transcript -> gene direct mapping
    for feat in features:
        if feat['type'] in transcript_types and feat['id'] and feat['parents']:
            for pid in feat['parents']:
                if pid in gene_ids:
                    transcript_to_gene[feat['id']] = pid
                    break

    # Second pass: transcript -> transcript -> gene chained mapping
    changed = True
    while changed:
        changed = False
        for feat in features:
            if feat['type'] not in transcript_types or not feat['id'] or feat['id'] in transcript_to_gene:
                continue
            for pid in feat['parents']:
                if pid in transcript_to_gene:
                    transcript_to_gene[feat['id']] = transcript_to_gene[pid]
                    changed = True
                    break

    # Assign each feature to a gene_id
    for feat in features:
        if feat['type'] == 'gene':
            feat['gene_id'] = feat['id']
            continue

        for pid in feat['parents']:
            if pid in gene_ids:
                feat['gene_id'] = pid
                break
            if pid in transcript_to_gene:
                feat['gene_id'] = transcript_to_gene[pid]
                break

    return features


def build_sequence_layout(dataframe, name_col, length_col):
    """Build cumulative coordinate layout for sequences"""
    names = list(dict.fromkeys(dataframe[name_col].astype(str).tolist()))
    length_map = dataframe.groupby(name_col)[length_col].max().to_dict()

    offsets = {}
    cursor = 0
    for name in names:
        offsets[name] = cursor
        cursor += int(length_map.get(name, 0))

    return names, length_map, offsets, cursor


def build_axis_layout(dataframe, name_col, length_col, region=None):
    """Build display layout; if region is set, only keep that local window."""
    normalized = normalize_region(region)
    if normalized is None:
        return build_sequence_layout(dataframe, name_col, length_col)

    chrom, start, end = normalized
    chrom_df = dataframe[dataframe[name_col].astype(str) == chrom]
    if chrom_df.empty:
        raise ValueError(f'Region chromosome {chrom} not found in {name_col}.')

    seq_len = int(chrom_df[length_col].max())
    clipped_start = max(0, min(start, seq_len))
    clipped_end = max(0, min(end, seq_len))
    if clipped_end <= clipped_start:
        raise ValueError(f'Region {chrom}:{start}-{end} is outside sequence range [0, {seq_len}].')

    local_len = clipped_end - clipped_start
    return [chrom], {chrom: local_len}, {chrom: 0}, local_len


def get_axis_scale(total_length_bp):
    """Determine axis scale based on total displayed length"""
    return 1e3 if total_length_bp < 1e6 else 1e6


def get_axis_label(axis_scale):
    """Get axis unit label"""
    return 'Mb' if axis_scale >= 1e6 else 'Kb'


def _points_to_figure_fraction(fig, points, axis='height'):
    """Convert point units to figure-relative width/height fraction."""
    fig_inches = fig.get_figheight() if axis == 'height' else fig.get_figwidth()
    if fig_inches <= 0:
        return 0.0
    return (float(points) / 72.0) / fig_inches


def map_interval_to_global(chrom, start, end, seq_offsets, region=None):
    """Map sequence-local interval to global concatenated coordinates"""
    if chrom not in seq_offsets:
        logging.warning('No chrom exsists in seqoffsets.')
        return None

    draw_start = float(start)
    draw_end = float(end)

    normalized = normalize_region(region)
    r_start = 0
    if normalized is not None:
        r_chrom, r_start, r_end = normalized
        if chrom != r_chrom:
            return None
        draw_start = max(draw_start, float(r_start))
        draw_end = min(draw_end, float(r_end))

    if draw_end <= draw_start:
        return None

    # If a local region is enabled, convert genomic coordinates into
    # region-local coordinates so the displayed axis starts from 0.
    if normalized is not None:
        draw_start -= float(r_start)
        draw_end -= float(r_start)

    return seq_offsets[chrom] + draw_start, seq_offsets[chrom] + draw_end


def _intersect_line_with_interval(v0, v1, low, high, u0, u1):
    """Intersect parametric line value v(u)=v0+(v1-v0)u with [low, high]."""
    dv = float(v1) - float(v0)
    low = float(low)
    high = float(high)

    if dv == 0:
        if low <= float(v0) <= high:
            return u0, u1
        return None

    c1 = (low - float(v0)) / dv
    c2 = (high - float(v0)) / dv
    lo = min(c1, c2)
    hi = max(c1, c2)

    u0_new = max(float(u0), lo)
    u1_new = min(float(u1), hi)
    if u1_new <= u0_new:
        return None
    return u0_new, u1_new


def clip_alignment_pair_by_regions(
    query_start,
    query_end,
    target_start,
    target_end,
    strand,
    query_region=None,
    target_region=None,
):
    """Coupled clipping for one alignment segment.

    Keep query/target endpoints synchronized in one shared parameter space,
    so clipping on either axis will update the opposite axis consistently.
    """
    qs = float(query_start)
    qe = float(query_end)
    ts = float(target_start)
    te = float(target_end)

    if qe <= qs or te <= ts:
        return None

    # Parametric segment in [0, 1] along query direction.
    u0, u1 = 0.0, 1.0

    # Clip by query region (if any).
    normalized_q = normalize_region(query_region) if query_region is not None else None
    if normalized_q is not None:
        _, q_low, q_high = normalized_q
        clipped = _intersect_line_with_interval(qs, qe, q_low, q_high, u0, u1)
        if clipped is None:
            return None
        u0, u1 = clipped

    # Build target line orientation used by plotting.
    if strand == '-':
        t0, t1 = te, ts
    else:
        t0, t1 = ts, te

    # Clip by target region (if any).
    normalized_t = normalize_region(target_region) if target_region is not None else None
    if normalized_t is not None:
        _, t_low, t_high = normalized_t
        clipped = _intersect_line_with_interval(t0, t1, t_low, t_high, u0, u1)
        if clipped is None:
            return None
        u0, u1 = clipped

    # Recover clipped query interval.
    q_clip_start = qs + (qe - qs) * u0
    q_clip_end = qs + (qe - qs) * u1

    # Recover clipped target interval from oriented target line,
    # then convert to ascending genomic coordinates.
    t_clip_a = t0 + (t1 - t0) * u0
    t_clip_b = t0 + (t1 - t0) * u1
    t_clip_start = min(t_clip_a, t_clip_b)
    t_clip_end = max(t_clip_a, t_clip_b)

    if q_clip_end <= q_clip_start or t_clip_end <= t_clip_start:
        return None

    return q_clip_start, q_clip_end, t_clip_start, t_clip_end


def draw_bed_rectangles(ax, bed_data, seq_offsets, orientation='x', region=None, axis_scale=1e6):
    """Draw background rectangles from BED file in global coordinates"""
    if not bed_data:
        return

    for chrom, start, end, color in bed_data:
        if chrom not in seq_offsets:
            continue
        mapped = map_interval_to_global(chrom, start, end, seq_offsets, region=region)
        if mapped is None:
            continue
        global_start, global_end = mapped

        rect_color = color if color else COLORS['rect_default']

        try:
            from matplotlib.colors import to_rgba
            rgba = to_rgba(rect_color)
            rect_color_with_alpha = (rgba[0], rgba[1], rgba[2], RECT_ALPHA)
        except Exception:
            rect_color_with_alpha = rect_color

        if orientation == 'x':
            ax.axvspan(global_start / axis_scale, global_end / axis_scale,
                       color=rect_color_with_alpha, zorder=0)
        else:
            ax.axhspan(global_start / axis_scale, global_end / axis_scale,
                       color=rect_color_with_alpha, zorder=0)


def draw_bed_markers(
    ax,
    bed_data,
    seq_offsets,
    orientation='x',
    region=None,
    axis_scale=1e6,
    gap_axes=0.0,
    tri_size_pt=TRIANGLE_MARKER_SIZE_PT,
    tri_offset_pt=0.0,
    marker_ax=None,
    draw_guides=True,
):
    """Draw triangle markers and guide dashed lines.

    Triangle size and offset use point units (via scatter + ScaledTranslation),
    so marker appearance is stable across global/local zoom ranges.
    """
    if not bed_data:
        return

    for chrom, start, end, color in bed_data:
        if chrom not in seq_offsets:
            continue
        mapped = map_interval_to_global(chrom, start, end, seq_offsets, region=region)
        if mapped is None:
            continue
        global_start, global_end = mapped
        midpoint = (global_start + global_end) / 2 / axis_scale
        tri_color = _safe_color(color, COLORS['triangle_default'])

        if orientation == 'x':
            # Dashed guide line in main alignment panel
            if draw_guides and ax is not None:
                ax.axvline(
                    midpoint,
                    color=tri_color,
                    linestyle=(0, (4, 3)),
                    linewidth=0.9,
                    alpha=0.8,
                    zorder=1,
                )

            if marker_ax is None:
                # x in data, y in axes fraction; marker size is in points^2.
                # Add a fixed point offset along y so distance to axis is zoom-independent.
                trans = ax.get_xaxis_transform() + mtransforms.ScaledTranslation(
                    0,
                    tri_offset_pt / 72.0,
                    ax.figure.dpi_scale_trans,
                )
                ax.scatter(
                    [midpoint],
                    [1 + gap_axes],
                    s=tri_size_pt ** 2,
                    marker='v',
                    c=[tri_color],
                    transform=trans,
                    clip_on=False,
                    zorder=10,
                )
            else:
                marker_ax.scatter(
                    [midpoint],
                    [0.5],
                    s=tri_size_pt ** 2,
                    marker='v',
                    c=[tri_color],
                    clip_on=False,
                    zorder=10,
                )
        else:
            # Dashed guide line in main alignment panel
            if draw_guides and ax is not None:
                ax.axhline(
                    midpoint,
                    color=tri_color,
                    linestyle=(0, (4, 3)),
                    linewidth=0.9,
                    alpha=0.8,
                    zorder=1,
                )

            if marker_ax is None:
                # x in axes fraction, y in data; marker size is in points^2.
                # Add a fixed point offset along x so distance to axis is zoom-independent.
                trans = ax.get_yaxis_transform() + mtransforms.ScaledTranslation(
                    tri_offset_pt / 72.0,
                    0,
                    ax.figure.dpi_scale_trans,
                )
                ax.scatter(
                    [1 + gap_axes],
                    [midpoint],
                    s=tri_size_pt ** 2,
                    marker='>',
                    c=[tri_color],
                    transform=trans,
                    clip_on=False,
                    zorder=10,
                )
            else:
                marker_ax.scatter(
                    [0.5],
                    [midpoint],
                    s=tri_size_pt ** 2,
                    marker='>',
                    c=[tri_color],
                    clip_on=False,
                    zorder=10,
                )


def _calculate_nice_tick_step(
    length_values_bp: List[int], 
    target_ticks_per_seq: int = 3
) -> int:
    """Find appropraite tick span."""
    valid_lengths = [int(v) for v in length_values_bp if int(v) > 0]
    if not valid_lengths:
        return int(1e6)
    
    avg_len_bp = sum(valid_lengths) / len(valid_lengths)
    raw_step_bp = max(1.0, avg_len_bp / max(1, int(target_ticks_per_seq)))
    
    exponent = math.floor(math.log10(raw_step_bp)) if raw_step_bp > 0 else 0
    candidates = []
    for exp in (exponent - 1, exponent, exponent + 1, exponent + 2):
        base = 10 ** exp
        for multiplier in (1, 2, 5):
            candidates.append(multiplier * base)

    best_step = min(candidates, key=lambda c: abs(c - raw_step_bp))
    return max(1, int(round(best_step)))


def _format_tick_label(units_value: float) -> str:
    if abs(units_value - round(units_value)) < 1e-9:
        return str(int(round(units_value)))
    return f'{units_value:.2f}'.rstrip('0').rstrip('.')


def _build_sequence_axis_ticks(
    seq_names: List[str],
    seq_lengths: Dict[str, int],
    seq_offsets: Dict[str, int],
    step_bp: int,
    scale_bp: float,
    axis_region: Optional[Tuple[str, int, int]] = None
) -> Tuple[List[float], List[str]]:
    """Build ticks and labels for axis."""
    ticks = []
    labels = []
    
    region_start_bp = None
    if axis_region is not None:
        normalized = normalize_region(axis_region)
        if normalized is not None:
            _, region_start_bp, _ = normalized
    
    for seq in seq_names:
        seq_len = int(seq_lengths.get(seq, 0))
        seq_off = int(seq_offsets.get(seq, 0))
        if seq_len <= 0:
            continue

        local_bp = 0
        while local_bp < seq_len:
            ticks.append((seq_off + local_bp) / scale_bp)
            
            if local_bp == 0:
                if region_start_bp is None:
                    labels.append(str(seq))
                else:
                    labels.append(_format_tick_label(region_start_bp / scale_bp))
            else:
                if region_start_bp is None:
                    labels.append(_format_tick_label(local_bp / scale_bp))
                else:
                    labels.append(_format_tick_label((region_start_bp + local_bp) / scale_bp))
            
            local_bp += step_bp

        if local_bp == 0:
            ticks.append(seq_off / scale_bp)
            labels.append('0')
    
    return ticks, labels


def _draw_sequence_boundary_lines(
    ax,
    x_names: List[str],
    x_offsets: Dict[str, int],
    y_names: List[str],
    y_offsets: Dict[str, int],
    axis_scale: float
) -> None:
    for x_name in x_names:
        start_bp = x_offsets[x_name]
        if start_bp > 0:
            ax.axvline(
                start_bp / axis_scale, 
                color='lightgray', 
                linestyle='--', 
                linewidth=0.8, 
                zorder=0
            )    

    for y_name in y_names:
        start_bp = y_offsets[y_name]
        if start_bp > 0:
            ax.axhline(
                start_bp / axis_scale, 
                color='lightgray', 
                linestyle='--', 
                linewidth=0.8, 
                zorder=0
            )


def _configure_axis_ticks_style(ax) -> None:
    """Set content for axis ticks"""
    ax.tick_params(
        axis='x',
        which='major',
        bottom=True,
        top=False,
        labelbottom=True,
        labeltop=False,
        length=2,
        labelsize=LABEL_SIZE,
    )
    ax.tick_params(
        axis='y',
        which='major',
        left=True,
        right=False,
        labelleft=True,
        labelright=False,
        length=2,
        labelsize=LABEL_SIZE,
    )
    
    # Rotate x label
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(90)
        lbl.set_fontsize(LABEL_SIZE)


def draw_sequence_boundaries_and_labels(
    ax,
    x_layout: AxisLayout,
    y_layout: AxisLayout,
    axis_scale: float,
) -> None:
    """Set boudaries and labels for sequences
    1. draw sequence boundaries
    2. set ticks for axis in that sequence tick starting from 0
    3. replace 0 with sequence name
    4. adjust the density of ticks depending on region restriction
    """
    x_names = x_layout.names
    x_lengths = x_layout.lengths
    x_offsets = x_layout.offsets
    x_region = x_layout.region

    y_names = y_layout.names
    y_lengths = y_layout.lengths
    y_offsets = y_layout.offsets
    y_region = y_layout.region

    # Draw sequence boundaries
    _draw_sequence_boundary_lines(
        ax, x_names, x_offsets, y_names, y_offsets, axis_scale
    )
    
    # Calculate tick span
    x_target_ticks = 5 if x_region is not None else 3
    y_target_ticks = 5 if y_region is not None else 3
    
    x_step_bp = _calculate_nice_tick_step(
        list(x_lengths.values()), 
        target_ticks_per_seq=x_target_ticks
    )
    y_step_bp = _calculate_nice_tick_step(
        list(y_lengths.values()), 
        target_ticks_per_seq=y_target_ticks
    )
    
    # Set ticks and labels
    x_ticks, x_labels = _build_sequence_axis_ticks(
        x_names,
        x_lengths,
        x_offsets,
        x_step_bp,
        axis_scale,
        axis_region=x_region,
    )
    y_ticks, y_labels = _build_sequence_axis_ticks(
        y_names,
        y_lengths,
        y_offsets,
        y_step_bp,
        axis_scale,
        axis_region=y_region,
    )
    
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    
    _configure_axis_ticks_style(ax)

def _has_mappable_bed_data(bed_data, seq_offsets, region=None):
    """Return True if at least one BED interval is drawable on current axis layout."""
    if not bed_data:
        return False

    for chrom, start, end, _ in bed_data:
        if chrom not in seq_offsets:
            continue
        mapped = map_interval_to_global(chrom, start, end, seq_offsets, region=region)
        if mapped is not None:
            return True
    return False


def _has_drawable_gene_data(genes, seq_offsets, region=None):
    """Return True if at least one gene (or pseudo-gene) is drawable."""
    if not genes:
        return False

    cds_list, exon_list, gene_list = _filter_gene_features(genes)
    gene_list = _merge_gene_features(cds_list, exon_list, gene_list)

    for gene in gene_list:
        chrom = gene.get('chrom')
        if chrom not in seq_offsets:
            continue
        mapped = map_interval_to_global(chrom, gene['start'], gene['end'], seq_offsets, region=region)
        if mapped is not None:
            return True

    return False


def estimate_gene_track_count(genes, seq_offsets, region=None, max_tracks=40):
    """Estimate how many non-overlapping gene tracks are needed."""
    if not genes:
        return 0

    cds_list, exon_list, gene_list = _filter_gene_features(genes)
    gene_list = _merge_gene_features(cds_list, exon_list, gene_list)

    drawable_intervals = []
    for gene in gene_list:
        chrom = gene.get('chrom')
        if chrom not in seq_offsets:
            continue
        mapped = map_interval_to_global(chrom, gene['start'], gene['end'], seq_offsets, region=region)
        if mapped is None:
            continue
        drawable_intervals.append(mapped)

    if not drawable_intervals:
        return 0

    drawable_intervals.sort(key=lambda x: (x[0], x[1]))
    track_ends = []
    for g_start, g_end in drawable_intervals:
        assigned = False
        for tid, last_end in enumerate(track_ends):
            if g_start > last_end:
                track_ends[tid] = g_end
                assigned = True
                break
        if not assigned:
            track_ends.append(g_end)

    return min(max(1, len(track_ends)), max_tracks)

def _merge_gene_features_from_subfeatures(cds_list, exon_list, gene_list):
    """从CDS/exon构建伪基因范围 (如果没有显式gene记录)
    
    Args:
        cds_list: CDS特征列表
        exon_list: exon特征列表
        gene_list: gene特征列表
        
    Returns:
        合并后的gene列表
    """
    if gene_list:
        return gene_list
    
    merged = {}
    for feat in cds_list + exon_list:
        gid = feat['gene_id']
        if not gid:
            continue
        if gid not in merged:
            merged[gid] = {
                'chrom': feat['chrom'],
                'type': 'gene',
                'start': feat['start'],
                'end': feat['end'],
                'strand': feat['strand'],
                'gene_id': gid,
            }
        else:
            merged[gid]['start'] = min(merged[gid]['start'], feat['start'])
            merged[gid]['end'] = max(merged[gid]['end'], feat['end'])
    
    return list(merged.values())


def _prepare_drawable_genes(gene_list, seq_offsets, region, axis_scale):
    """准备可绘制的基因数据
    
    Args:
        gene_list: 基因特征列表
        seq_offsets: 序列偏移字典
        region: 区域限制
        axis_scale: 轴刻度
        
    Returns:
        可绘制基因列表
    """
    draw_genes = []
    for gene in gene_list:
        gid = gene.get('gene_id') or gene.get('id')
        if not gid:
            continue
        if gene['chrom'] not in seq_offsets:
            continue
        
        mapped = map_interval_to_global(
            gene['chrom'], 
            gene['start'], 
            gene['end'], 
            seq_offsets, 
            region=region
        )
        if mapped is None:
            continue
        
        g_start, g_end = mapped
        draw_genes.append({
            'gene_id': gid,
            'chrom': gene['chrom'],
            'strand': gene.get('strand', '.'),
            'global_start': g_start,
            'global_end': g_end,
            'color': gene.get('color'),
            'show_gene_id': bool(gene.get('show_gene_id', False)),
        })
    
    return draw_genes


def _assign_genes_to_tracks(draw_genes, max_tracks=40):
    """Allocate non-overlap tracks for genes
    Returns:
        (gene lists, true tracks)
    """
    if not draw_genes:
        return [], 0

    draw_genes.sort(key=lambda x: (x['global_start'], x['global_end']))

    track_ends = []
    for gene in draw_genes:
        assigned = False
        for tid, last_end in enumerate(track_ends):
            if gene['global_start'] > last_end:
                gene['track'] = tid
                track_ends[tid] = gene['global_end']
                assigned = True
                break
        if not assigned:
            gene['track'] = len(track_ends)
            track_ends.append(gene['global_end'])
    
    n_tracks = min(max(1, len(track_ends)), max_tracks)
    return draw_genes, n_tracks


def _setup_gene_annotation_axis(ax, orientation, n_tracks, draw_ax=None):
    """设置基因注释轴
    
    Args:
        ax: 主轴对象
        orientation: 方向 ('top' 或 'right')
        n_tracks: 轨道数量
        draw_ax: 预定义的绘图轴 (可选)
        
    Returns:
        注释轴对象
    """
    if draw_ax is not None:
        return draw_ax
    
    ax_pos = ax.get_position()
    if orientation == 'top':
        inset_ax = ax.figure.add_axes([ax_pos.x0, ax_pos.y1, ax_pos.width, 0.08])
    else:  # 'right'
        inset_ax = ax.figure.add_axes([ax_pos.x1, ax_pos.y0, 0.08, ax_pos.height])
    
    return inset_ax


def _configure_gene_axis_limits(inset_ax, ax, orientation, n_tracks):
    track_margin = max((GENE_TRACK_PITCH - LAYOUT['max_feat_span']) / 2, LAYOUT['track_margin_min'])
    max_track_pos = (n_tracks - 1) * GENE_TRACK_PITCH
    
    if orientation == 'top':
        inset_ax.set_xlim(ax.get_xlim())
        inset_ax.set_ylim(-track_margin, max_track_pos + track_margin)
    else:  # 'right'
        inset_ax.set_ylim(ax.get_ylim())
        inset_ax.set_xlim(-track_margin, max_track_pos + track_margin)

    inset_ax.set_xticks([])
    inset_ax.set_yticks([])
    for spine in inset_ax.spines.values():
        spine.set_visible(False)
    inset_ax.set_facecolor('none')


def _get_subfeatures_for_gene(gene_id, cds_by_gene, exon_by_gene):
    if gene_id in cds_by_gene:
        return cds_by_gene[gene_id], 0.55, COLORS['gene_cds']
    elif gene_id in exon_by_gene:
        return exon_by_gene[gene_id], 0.35, COLORS['gene_exon']
    else:
        return [], 0.35, COLORS['gene_exon']


def _draw_horizontal_gene_track(
    inset_ax, gene, cds_by_gene, exon_by_gene, 
    seq_offsets, region, axis_scale, ax
):
    track_pos = gene['track'] * GENE_TRACK_PITCH
    g1 = gene['global_start'] / axis_scale
    g2 = gene['global_end'] / axis_scale
    strand = gene['strand']
    gid = gene['gene_id']
    gene_color = _safe_color(gene.get('color'), 'black')
    
    # Draw gene backbones
    inset_ax.plot(
        [g1, g2], 
        [track_pos, track_pos], 
        '-', 
        linewidth=1.0, 
        color=gene_color, 
        zorder=2
    )
    
    # Draw directions arrow of genes
    x_delta = max((ax.get_xlim()[1] - ax.get_xlim()[0]) * LAYOUT['delta_factor'], 1e-6)
    if strand == '+':
        inset_ax.annotate(
            '', 
            xy=(g2, track_pos), 
            xytext=(max(g1, g2 - x_delta), track_pos),
            arrowprops=dict(arrowstyle='->', color=gene_color, lw=1)
        )
    elif strand == '-':
        inset_ax.annotate(
            '', 
            xy=(g1, track_pos), 
            xytext=(min(g2, g1 + x_delta), track_pos),
            arrowprops=dict(arrowstyle='->', color=gene_color, lw=1)
        )
    
    # Draw sub-features
    sub_feats, feat_height, default_color = _get_subfeatures_for_gene(
        gid, cds_by_gene, exon_by_gene
    )
    
    for feat in sub_feats:
        if feat['chrom'] != gene['chrom']:
            continue
        
        mapped = map_interval_to_global(
            feat['chrom'], 
            feat['start'], 
            feat['end'], 
            seq_offsets, 
            region=region
        )
        if mapped is None:
            continue
        
        f1, f2 = mapped
        rect = mpatches.Rectangle(
            (f1 / axis_scale, track_pos - feat_height / 2),
            max((f2 - f1) / axis_scale, 1e-9),
            feat_height,
            color=_safe_color(feat.get('color'), default_color),
            zorder=3,
        )
        inset_ax.add_patch(rect)
    
    # Draw gene labels
    if gene.get('show_gene_id'):
        x_label_pad = max((ax.get_xlim()[1] - ax.get_xlim()[0]) * LAYOUT['label_pad_factor'], 1e-6)
        inset_ax.text(
            g2 + x_label_pad,
            track_pos,
            gid,
            fontsize=GENE_SIZE,
            va='center',
            ha='left',
            color=gene_color,
            clip_on=False,
            zorder=5,
        )


def _draw_vertical_gene_track(
    inset_ax, gene, cds_by_gene, exon_by_gene, 
    seq_offsets, region, axis_scale, ax
):
    track_pos = gene['track'] * GENE_TRACK_PITCH
    g1 = gene['global_start'] / axis_scale
    g2 = gene['global_end'] / axis_scale
    strand = gene['strand']
    gid = gene['gene_id']
    gene_color = _safe_color(gene.get('color'), 'black')
    
    # Draw backbone of genes
    inset_ax.plot(
        [track_pos, track_pos], 
        [g1, g2], 
        '-', 
        linewidth=1.0, 
        color=gene_color, 
        zorder=2
    )
    
    # Draw arrow of direction
    y_delta = max((ax.get_ylim()[1] - ax.get_ylim()[0]) * LAYOUT['delta_factor'], 1e-6)
    if strand == '+':
        inset_ax.annotate(
            '', 
            xy=(track_pos, g2), 
            xytext=(track_pos, max(g1, g2 - y_delta)),
            arrowprops=dict(arrowstyle='->', color=gene_color, lw=1)
        )
    elif strand == '-':
        inset_ax.annotate(
            '', 
            xy=(track_pos, g1), 
            xytext=(track_pos, min(g2, g1 + y_delta)),
            arrowprops=dict(arrowstyle='->', color=gene_color, lw=1)
        )
    
    # Draw sub-features
    sub_feats, feat_width, default_color = _get_subfeatures_for_gene(
        gid, cds_by_gene, exon_by_gene
    )
    
    for feat in sub_feats:
        if feat['chrom'] != gene['chrom']:
            continue
        
        mapped = map_interval_to_global(
            feat['chrom'], 
            feat['start'], 
            feat['end'], 
            seq_offsets, 
            region=region
        )
        if mapped is None:
            continue
        
        f1, f2 = mapped
        rect = mpatches.Rectangle(
            (track_pos - feat_width / 2, f1 / axis_scale),
            feat_width,
            max((f2 - f1) / axis_scale, 1e-9),
            color=_safe_color(feat.get('color'), default_color),
            zorder=3,
        )
        inset_ax.add_patch(rect)
    
    # 绘制基因ID标签
    if gene.get('show_gene_id'):
        track_label_pad = LAYOUT['track_label_pad']
        inset_ax.text(
            track_pos + track_label_pad,
            (g1 + g2) / 2,
            gid,
            fontsize=GENE_SIZE,
            va='center',
            ha='left',
            color=gene_color,
            clip_on=False,
            zorder=5,
        )


def draw_gene_structure(
    ax,
    genes,
    context: GeneDrawContext,
    orientation='top',
    max_tracks=40,
    draw_ax=None,
):
    """Draw track of CDS/exon"""

    seq_offsets = context.seq_offsets
    axis_scale = context.axis_scale
    region = context.region

    if not genes:
        return False

    cds_list = [g for g in genes if g['type'] == 'CDS' and g['gene_id']]
    exon_list = [g for g in genes if g['type'] == 'exon' and g['gene_id']]
    gene_list = [g for g in genes if g['type'] == 'gene']

    gene_list = _merge_gene_features_from_subfeatures(cds_list, exon_list, gene_list)

    cds_by_gene = {}
    for cds in cds_list:
        cds_by_gene.setdefault(cds['gene_id'], []).append(cds)
    
    exon_by_gene = {}
    for exon in exon_list:
        exon_by_gene.setdefault(exon['gene_id'], []).append(exon)

    draw_genes = _prepare_drawable_genes(gene_list, seq_offsets, region, axis_scale)
    if not draw_genes:
        return False

    draw_genes, n_tracks = _assign_genes_to_tracks(draw_genes, max_tracks)

    inset_ax = _setup_gene_annotation_axis(ax, orientation, n_tracks, draw_ax)
    _configure_gene_axis_limits(inset_ax, ax, orientation, n_tracks)
    
    # Draw all genes
    for gene in draw_genes:
        if gene['track'] >= n_tracks:
            continue
        
        if orientation == 'top':
            _draw_horizontal_gene_track(
                inset_ax, gene, cds_by_gene, exon_by_gene,
                seq_offsets, region, axis_scale, ax
            )
        else:  # 'right'
            _draw_vertical_gene_track(
                inset_ax, gene, cds_by_gene, exon_by_gene,
                seq_offsets, region, axis_scale, ax
            )
    
    return True


def draw_dotplot_with_highlight(
    ctx: PlotRuntimeContext,
):
    """Draw one global dotplot including all query-target sequences"""
    dataframe = ctx.dataframe
    if dataframe.empty:
        logging.warning('No alignments to draw after filtering.')
        return

    query_layout = ctx.query_layout
    target_layout = ctx.target_layout
    query_region = query_layout.region
    target_region = target_layout.region

    query_offsets = query_layout.offsets
    target_offsets = target_layout.offsets
    query_total = query_layout.total
    target_total = target_layout.total

    query_genes = ctx.query_genes
    target_genes = ctx.target_genes
    beds = ctx.beds

    if target_total <= 0 or query_total <= 0:
        logging.warning('Invalid target/query total length; skip drawing.')
        return

    axis_scale = ctx.axis_scale
    unit_label = ctx.unit_label

    has_top_tri = _has_mappable_bed_data(beds.tri_x, query_offsets, region=query_region)
    has_right_tri = _has_mappable_bed_data(beds.tri_y, target_offsets, region=target_region)
    has_top_gene = _has_drawable_gene_data(query_genes, query_offsets, region=query_region)
    has_right_gene = _has_drawable_gene_data(target_genes, target_offsets, region=target_region)
    top_gene_tracks = estimate_gene_track_count(query_genes, query_offsets, region=query_region) if has_top_gene else 0
    right_gene_tracks = estimate_gene_track_count(target_genes, target_offsets, region=target_region) if has_right_gene else 0

    # Annotation track order: triangle first, then gene structure.
    top_layers = []
    if has_top_tri:
        top_layers.append('triangle')
    if has_top_gene:
        top_layers.append('gene')

    right_layers = []
    if has_right_tri:
        right_layers.append('triangle')
    if has_right_gene:
        right_layers.append('gene')

    fig, ax = plt.subplots(figsize=(8, 8))

    triangle_slot_h = max(
        _points_to_figure_fraction(
            fig,
            TRIANGLE_MARKER_SIZE_PT + 2 * TRIANGLE_TRACK_PAD_PT,
            axis='height',
        ),
        0.004,
    )
    triangle_slot_w = max(
        _points_to_figure_fraction(
            fig,
            TRIANGLE_MARKER_SIZE_PT + 2 * TRIANGLE_TRACK_PAD_PT,
            axis='width',
        ),
        0.004,
    )
    top_gene_slot_h = max(
        _points_to_figure_fraction(
            fig,
            max(1, top_gene_tracks) * GENE_TRACK_SLOT_PT + 2 * GENE_TRACK_PAD_PT,
            axis='height',
        ),
        0.008,
    )
    right_gene_slot_w = max(
        _points_to_figure_fraction(
            fig,
            max(1, right_gene_tracks) * GENE_TRACK_SLOT_PT + 2 * GENE_TRACK_PAD_PT,
            axis='width',
        ),
        0.008,
    )

    top_slot_h_map = {
        'triangle': triangle_slot_h,
        'gene': top_gene_slot_h,
    }
    right_slot_w_map = {
        'triangle': triangle_slot_w,
        'gene': right_gene_slot_w,
    }
    top_track_gap = 0.0015
    right_track_gap = 0.0015

    top_total = sum(top_slot_h_map[layer] for layer in top_layers) + max(0, len(top_layers) - 1) * top_track_gap
    right_total = sum(right_slot_w_map[layer] for layer in right_layers) + max(0, len(right_layers) - 1) * right_track_gap

    main_left = 0.08
    main_bottom = 0.08
    main_right = 0.95 - right_total
    main_top = 0.95 - top_total
    ax.set_position((main_left, main_bottom, main_right - main_left, main_top - main_bottom))

    # X-axis = Query, Y-axis = Target.
    # When region is enabled, coordinates are already mapped to local region.
    x_min, x_max = 0, query_total / axis_scale
    y_min, y_max = 0, target_total / axis_scale

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    top_axes = {}
    top_cursor = main_top
    for idx, layer in enumerate(top_layers):
        h = top_slot_h_map[layer]
        y0 = top_cursor
        layer_ax = fig.add_axes((main_left, y0, main_right - main_left, h))
        layer_ax.set_facecolor('none')
        layer_ax.set_xticks([])
        layer_ax.set_yticks([])
        for spine in layer_ax.spines.values():
            spine.set_visible(False)
        layer_ax.set_xlim(x_min, x_max)
        if layer == 'triangle':
            layer_ax.set_ylim(0, 1)
        top_axes[layer] = layer_ax
        top_cursor += h
        if idx < len(top_layers) - 1:
            top_cursor += top_track_gap

    right_axes = {}
    right_cursor = main_right
    for idx, layer in enumerate(right_layers):
        w = right_slot_w_map[layer]
        x0 = right_cursor
        layer_ax = fig.add_axes((x0, main_bottom, w, main_top - main_bottom))
        layer_ax.set_facecolor('none')
        layer_ax.set_xticks([])
        layer_ax.set_yticks([])
        for spine in layer_ax.spines.values():
            spine.set_visible(False)
        layer_ax.set_ylim(y_min, y_max)
        if layer == 'triangle':
            layer_ax.set_xlim(0, 1)
        right_axes[layer] = layer_ax
        right_cursor += w
        if idx < len(right_layers) - 1:
            right_cursor += right_track_gap

    # BED rectangles as background
    draw_bed_rectangles(ax, beds.rect_x, query_offsets, orientation='x', region=query_region, axis_scale=axis_scale)
    draw_bed_rectangles(ax, beds.rect_y, target_offsets, orientation='y', region=target_region, axis_scale=axis_scale)

    # Clip path for all alignment lines
    clip_rect = Rectangle(
        (x_min, y_min),
        x_max - x_min,
        y_max - y_min,
        transform=ax.transData,
        fill=False,
        edgecolor='black',
        linewidth=1,
        zorder=0,
    )
    ax.add_patch(clip_rect)

    # Draw alignments
    for row in dataframe.itertuples(index=False):
        if row.query_name not in query_offsets or row.target_name not in target_offsets:
            continue

        clipped_pair = clip_alignment_pair_by_regions(
            row.query_start,
            row.query_end,
            row.target_start,
            row.target_end,
            row.strand,
            query_region=query_region,
            target_region=target_region,
        )
        if clipped_pair is None:
            continue

        q_clip_start, q_clip_end, t_clip_start, t_clip_end = clipped_pair

        mapped_query = map_interval_to_global(
            row.query_name,
            q_clip_start,
            q_clip_end,
            query_offsets,
            region=query_region,
        )
        if mapped_query is None:
            continue

        mapped_target = map_interval_to_global(
            row.target_name,
            t_clip_start,
            t_clip_end,
            target_offsets,
            region=target_region,
        )
        if mapped_target is None:
            continue

        q1, q2 = mapped_query
        t1, t2 = mapped_target

        qs, qe = q1 / axis_scale, q2 / axis_scale
        ts, te = t1 / axis_scale, t2 / axis_scale

        if row.strand == '+':
            xs, ys = [qs, qe], [ts, te]
        else:
            xs, ys = [qs, qe], [te, ts]

        line_artist = ax.plot(
            xs,
            ys,
            linestyle='-',
            linewidth=1.0,
            color='black',
            alpha=0.85,
            solid_capstyle='round',
            zorder=2,
        )[0]
        line_artist.set_clip_path(clip_rect)

    # Sequence boundaries and labels
    draw_sequence_boundaries_and_labels(
        ax,
        query_layout,
        target_layout,
        axis_scale,
    )

    # BED triangles on top (fixed distance to axis)
    if has_top_tri:
        draw_bed_markers(
            ax,
            beds.tri_x,
            query_offsets,
            orientation='x',
            region=query_region,
            axis_scale=axis_scale,
            marker_ax=top_axes.get('triangle'),
            draw_guides=True,
        )
    if has_right_tri:
        draw_bed_markers(
            ax,
            beds.tri_y,
            target_offsets,
            orientation='y',
            region=target_region,
            axis_scale=axis_scale,
            marker_ax=right_axes.get('triangle'),
            draw_guides=True,
        )

    # GFF annotations
    if has_top_gene:
        query_gene_ctx = GeneDrawContext(
            seq_offsets=query_offsets,
            axis_scale=axis_scale,
            region=query_region,
        )
        draw_gene_structure(
            ax,
            query_genes,
            context=query_gene_ctx,
            orientation='top',
            draw_ax=top_axes.get('gene'),
        )

    if has_right_gene:
        target_gene_ctx = GeneDrawContext(
            seq_offsets=target_offsets,
            axis_scale=axis_scale,
            region=target_region,
        )
        draw_gene_structure(
            ax,
            target_genes,
            context=target_gene_ctx,
            orientation='right',
            draw_ax=right_axes.get('gene'),
        )

    ax.set_xlabel(f'Query ({unit_label})', fontsize=TITLE_SIZE)
    ax.set_ylabel(f'Target ({unit_label})', fontsize=TITLE_SIZE)
    ax.grid(False)
    # ax.set_aspect('equal')

    out_pdf = f'{ctx.output_prefix}.pdf'
    # out_png = f'{output_prefix}.png'
    plt.savefig(out_pdf, format='pdf', bbox_inches='tight')
    # plt.savefig(out_png, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    logging.info('Saved plot: %s', out_pdf)


def run(config: PlotInputConfig):
    """Main workflow using dataclass config/context."""

    beds = BedCollections.from_files(
        config.bed_tri_x,
        config.bed_tri_y,
        config.bed_rect_x,
        config.bed_rect_y,
    )

    column_names = [
        'query_name', 'query_length', 'query_start', 'query_end', 'strand',
        'target_name', 'target_length', 'target_start', 'target_end',
        'residue_matches', 'alignment_block_length', 'mapping_quality',
    ]

    df = pd.read_table(config.paf_file, header=None, names=column_names, usecols=range(12))
    df = df[df['alignment_block_length'] > config.filter_len].copy()

    if df.empty:
        logging.warning('No alignment records remain after length filtering: > %s', config.filter_len)
        return

    # Validate/normalize regions against PAF content so invalid requests fail fast.
    target_region = config.target_region
    query_region = config.query_region
    if target_region is not None:
        target_region = ensure_region_exists_in_paf(df, target_region, axis='target')
    if query_region is not None:
        query_region = ensure_region_exists_in_paf(df, query_region, axis='query')

    query_layout = AxisLayout.build(
        df,
        'query_name',
        'query_length',
        region=query_region,
    )
    target_layout = AxisLayout.build(
        df,
        'target_name',
        'target_length',
        region=target_region,
    )

    axis_scale = get_axis_scale(max(target_layout.total, query_layout.total))
    unit_label = get_axis_label(axis_scale)

    runtime_ctx = PlotRuntimeContext(
        dataframe=df,
        target_layout=target_layout,
        query_layout=query_layout,
        target_genes=parse_gff3(config.target_gff) if config.target_gff else [],
        query_genes=parse_gff3(config.query_gff) if config.query_gff else [],
        beds=beds,
        axis_scale=axis_scale,
        unit_label=unit_label,
        output_prefix=config.output_prefix,
    )

    logging.info('Start plotting global concatenated dotplot with %d alignments', len(df))
    draw_dotplot_with_highlight(runtime_ctx)


def main(paf_file,  target_region, query_region, target_gff, query_gff,
                    bed_tri_x, bed_tri_y, bed_rect_x, bed_rect_y,
                    filter_len, output_prefix, ):

    config = PlotInputConfig(
        paf_file=paf_file,
        target_region=target_region,
        query_region=query_region,
        target_gff=target_gff,
        query_gff=query_gff,
        bed_tri_x=bed_tri_x,
        bed_tri_y=bed_tri_y,
        bed_rect_x=bed_rect_x,
        bed_rect_y=bed_rect_y,
        filter_len=filter_len,
        output_prefix=output_prefix,
    )
    run(config)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Enhanced PAF global dotplot visualization with GFF annotation'
    )

    parser.add_argument('-i', '--paf', dest='paf', type=str, required=True,
                        help='Path of PAF file')

    parser.add_argument('-tr', '--target-region', dest='target_region', type=str, default=None,
                        help='Target region (format: chrom:start-end)')
    parser.add_argument('-qr', '--query-region', dest='query_region', type=str, default=None,
                        help='Query region (format: chrom:start-end)')

    parser.add_argument('-tg', '--target-gff', dest='target_gff', type=str, default=None,
                        help='Target sequence GFF3 file')
    parser.add_argument('-qg', '--query-gff', dest='query_gff', type=str, default=None,
                        help='Query sequence GFF3 file')

    parser.add_argument('-btx', '--bed-tri-x', dest='bed_tri_x', type=str, default=None,
                        help='BED file for X-axis triangle markers')
    parser.add_argument('-bty', '--bed-tri-y', dest='bed_tri_y', type=str, default=None,
                        help='BED file for Y-axis triangle markers')

    parser.add_argument('-brx', '--bed-rect-x', dest='bed_rect_x', type=str, default=None,
                        help='BED file for X-axis rectangle backgrounds')
    parser.add_argument('-bry', '--bed-rect-y', dest='bed_rect_y', type=str, default=None,
                        help='BED file for Y-axis rectangle backgrounds')

    parser.add_argument('-flen', '--filter-len', dest='flen', default=2000, type=int,
                        help='Filter alignments shorter than this length')
    parser.add_argument('-o', '--outpre', dest='outpre', default='all_query_target_dotplot', type=str,
                        help='Output plot files')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # target_region = parse_region(args.target_region)
    # query_region = parse_region(args.query_region)

    main(
        args.paf,
        # target_region,
        # query_region,
        args.target_region,
        args.query_region,
        args.target_gff,
        args.query_gff,
        args.bed_tri_x,
        args.bed_tri_y,
        args.bed_rect_x,
        args.bed_rect_y,
        args.flen,
        args.outpre,
    )

    # tets_paf = r'/mnt/e/Bio_analysis/PWS_lr/W2_hap1_CHM13v2m_chr15.paf'
    # bed_trix = r'/mnt/e/Bio_analysis/PWS_lr/W2_hap1_chr.gap.bed'
    # target_gff = r'/mnt/e/Bio_analysis/PWS_lr/test_15q1113.gff3'
    # target_region = 'chr15:17186630-31196487'
    # main(tets_paf, target_region,  None, target_gff, None, bed_trix, None, None, None, 2000, 'all_query_target_dotplot')
    # main(tets_paf, None,  None, target_gff, None, bed_trix, None, None, None, 2000, 'all_query_target_dotplot')