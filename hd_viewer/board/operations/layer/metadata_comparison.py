"""Descriptive comparison of observation metadata; expression counts are not used."""
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class MetadataComparison:
    summary: pd.DataFrame
    columns: tuple
    groups: tuple
    statistic: str


def numeric_columns(frame):
    return [str(name) for name in frame.columns
            if pd.api.types.is_numeric_dtype(frame[name].dtype)
            and not pd.api.types.is_bool_dtype(frame[name].dtype)
            and not pd.api.types.is_complex_dtype(frame[name].dtype)]


def metadata_presets(adata):
    """Column sets of the saved deconvolution runs: '<run>:fraction' / '<run>:abundance' ->
    dict(title, labels (column -> cell type), unit, fractions)."""
    records = dict(adata.uns.get('deconvolution_runs', {}))
    if not records and 'deconvolution' in adata.uns:
        records['latest'] = adata.uns['deconvolution']
    specs = {}
    numeric = set(numeric_columns(adata.obs))
    ordered = sorted(records.items(), key=lambda item: (1, int(item[0])) if str(item[0]).isdigit() else (0, str(item[0])))
    for run, record in ordered:
        for kind, title, unit in (('fraction', 'fractions', 'proportion'),
                                  ('abundance', 'abundance', str(record.get('abundance_units', 'abundance')))):
            pairs = [(str(col), str(label)) for col, label in zip(record.get(kind + '_columns', []), record.get('cell_types', []))
                     if col in numeric]
            if pairs:
                specs[f'{run}:{kind}'] = dict(title=f"Run {run}: {record.get('method', 'deconvolution')} {title}",
                                              labels=dict(pairs), unit=unit, fractions=kind == 'fraction')
    return specs


def compare_metadata(frame, masks, columns, statistic='mean'):
    """Summarize each column within each mask, retaining overlaps in every group.

    Non-finite values are omitted per column; n records the denominator used.
    Fractions are neither summed nor renormalized over the selected columns.
    """
    if statistic not in ('mean', 'median'):
        raise ValueError('Choose mean or median.')
    columns = list(columns)
    if len(masks) < 2:
        raise ValueError('Select at least two groups or clusters.')
    if not columns or len(set(columns)) != len(columns):
        raise ValueError('Select one or more distinct numeric metadata columns.')
    if set(columns) - set(numeric_columns(frame)):
        raise ValueError('Selected metadata columns must exist and contain numeric values.')
    values = frame[columns].to_numpy(dtype=float, na_value=np.nan)
    rows = []
    reducer = np.mean if statistic == 'mean' else np.median
    for name, mask in masks.items():
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (len(frame),):
            raise ValueError('Group masks must align with metadata observations.')
        total = int(mask.sum())
        if not total:
            raise ValueError(f'Group {name!r} contains no observations.')
        group = values[mask]
        for position, column in enumerate(columns):
            valid = group[:, position][np.isfinite(group[:, position])]
            rows.append(dict(column=column, group=name, value=float(reducer(valid)) if len(valid) else np.nan,
                             n=int(len(valid)), total=total))
    summary = pd.DataFrame(rows)
    if not summary.n.any():
        raise ValueError('No finite metadata values exist in the selected groups.')
    return MetadataComparison(summary, tuple(columns), tuple(masks), statistic)


def rest_mask(defined_masks, n_observations):
    """Native observations outside every defined custom group, irrespective of visibility."""
    assigned = np.zeros(n_observations, dtype=bool)
    for mask in defined_masks.values():
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != assigned.shape:
            raise ValueError('Defined group masks must align with metadata observations.')
        assigned |= mask
    return ~assigned


def compare_metadata_to_rest(frame, masks, defined_masks, columns, statistic='mean'):
    """One comparison per selected group/cluster against its own remainder.

    The rest of a selected set is every observation outside all defined custom groups and outside
    the set itself: as for DE vs Rest, a set never shares observations with its rest (a cluster
    overlapping no group would otherwise be part of its own comparison baseline).
    """
    if not masks:
        raise ValueError('Select at least one group or cluster to compare with The rest.')
    remainder = rest_mask(defined_masks, len(frame))
    if not remainder.any():
        raise ValueError('The rest contains no observations outside the defined groups.')
    rest_name = 'The rest'
    if rest_name in masks or rest_name in defined_masks:
        rest_name = 'The rest (outside defined groups)'
        while rest_name in masks or rest_name in defined_masks:
            rest_name += ' (remainder)'
    results = {}
    for name, mask in masks.items():
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != remainder.shape:
            raise ValueError('Group masks must align with metadata observations.')
        rest = remainder & ~mask
        if not rest.any():
            raise ValueError(f'The rest contains no observations outside the defined groups and {name!r}.')
        results[name] = compare_metadata(frame, {name: mask, rest_name: rest}, columns, statistic)
    return results


def build_metadata_plot(result, labels=None, colors=None, unit='value', fractions=False):
    """Grouped horizontal bars, with one bar per group for each metadata field."""
    from bokeh.models import ColumnDataSource, HoverTool, Range1d
    from bokeh.plotting import figure
    from bokeh.transform import dodge
    from ...utils.palette import GROUP_PALETTE, palette_color

    labels, colors = labels or {}, colors or {}
    names = [labels.get(column, column) for column in result.columns]
    if len(set(names)) != len(names):
        names = list(result.columns)
    mapping = dict(zip(result.columns, names))
    finite = result.summary.value.to_numpy()
    finite = finite[np.isfinite(finite)]
    low, high = min(0., float(finite.min())), max(0., float(finite.max()))
    margin = .05 * (high - low or 1.)
    x_range = Range1d(0, 1) if fractions else Range1d(low - (margin if low < 0 else 0), high + margin)
    chart = figure(y_range=list(reversed(names)), x_range=x_range,
                   height=max(320, min(2200, 34 * len(names) + 100)), sizing_mode='stretch_width',
                   title='Metadata comparison', tools='pan,wheel_zoom,box_zoom,reset,save')
    chart.xaxis.axis_label = f'{result.statistic.title()} {unit}'
    chart.yaxis.axis_label = 'Cell type' if fractions else 'Metadata column'
    width = .8 / len(result.groups)
    for position, name in enumerate(result.groups):
        frame = result.summary[result.summary.group == name].copy()
        frame['label'] = frame.column.map(mapping)
        frame['left'] = np.minimum(frame.value, 0)
        frame['right'] = np.maximum(frame.value, 0)
        source = ColumnDataSource(frame)
        chart.hbar(y=dodge('label', -.4 + width / 2 + position * width, range=chart.y_range),
                   left='left', right='right', height=width * .9, source=source,
                   color=colors.get(name) or palette_color(GROUP_PALETTE, position), legend_label=name)
    chart.add_tools(HoverTool(tooltips=[('Column', '@label'), ('Group', '@group'),
                                       (result.statistic.title(), '@value{0.0000}'),
                                       ('Finite observations', '@n'), ('Group observations', '@total')]))
    chart.legend.click_policy = 'hide'
    chart.legend.location = 'top_right'
    chart.ygrid.grid_line_color = None
    return chart
