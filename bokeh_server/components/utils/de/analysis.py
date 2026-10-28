
from .... import settings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from src.utils.progress import report

@dataclass
class DEResult:
    """Statistical results and the aligned expression/membership data used in testing.

    Expression is log1p-transformed after filtering genes. Its matrix is shared
    with the analysis cache and should be treated as read-only by consumers.
    """

    results: dict[str, pd.DataFrame]
    expression: ad.AnnData
    masks: dict[str, np.ndarray]
    names: tuple[str, ...]
    mode: int
    source: str
    factor: int
    n_overlap: int


class DEAnalysis:
    """Differential expression over explicit named observation masks.

    Instantiate with a sample exposing ``adata``; aggregated analysis additionally
    uses the sample's bin aggregation API. ``run`` is synchronous and can be used
    without GroupsManager, a component, a Bokeh document, or output files.
    """

    def __init__(self, sample):
        self.sample = sample
        self.prepared_data = {}

    def run(self, masks, mode=0, *, source="groups", factor=1):
        """Return a DEResult for each named set versus rest (0) or selected sets (1).

        Masks align with the sample's native observations. With factor > 1, counts
        are summed into bin aggregates and membership requires > 50% of their bins.
        Inputs are copied so preprocessing never modifies the sample's counts.
        """
        if mode not in (0, 1):
            raise ValueError("DE mode must be 0 (vs rest) or 1 (between selected).")
        if factor < 1 or int(factor) != factor:
            raise ValueError("Aggregation factor must be a positive integer.")
        if not masks:
            raise ValueError("Select at least one set to compare.")
        masks = {name: np.asarray(mask, dtype=bool).copy() for name, mask in masks.items()}
        if any(mask.shape != (self.sample.adata.n_obs,) for mask in masks.values()):
            raise ValueError("Each DE mask must have one value per sample observation.")
        active = list(masks)
        factor = int(factor)
        unit, units = {"groups": ("group", "groups"), "clusters": ("cluster", "clusters")}.get(
            source, ("set", "groups/clusters"))
        obs_word = "observations" if factor <= 1 else f"{factor} × {factor} bins"
        if factor > 1:
            level = self.sample._get_aggregated_sample(factor).adata
            level_masks = self.masks_at_level(masks, factor, level.obs_names)
        else:
            level, level_masks = self.sample.adata, masks

        cache = self.prepared_data
        if factor not in cache:
            report(f"Preparing expression data ({level.n_obs} {obs_word})")
            _a = level.copy()
            sc.pp.filter_genes(_a, min_cells=settings.MIN_GENE_CELLS)
            sc.pp.log1p(_a)
            cache[factor] = _a

        # filter_genes/log1p keep the observations: masks still line up
        adata_copy = cache[factor]
        results = {}
        tested_data = adata_copy
        n_overlap = 0

        if mode == 0:
            for cname in active:
                if cname not in level_masks:
                    continue
                in_cluster = pd.Series(level_masks[cname], index=adata_copy.obs.index)
                if (~in_cluster).sum() == 0 or in_cluster.sum() < settings.DE_MIN_OBSERVATIONS:
                    report(f"Skipping {cname}: fewer than {settings.DE_MIN_OBSERVATIONS} {obs_word}")
                    continue
                report(f"Testing {cname} vs rest")
                adata_copy.obs['_de_group'] = in_cluster.map(
                    {True: cname, False: 'Rest'}
                ).astype(str)
                try:
                    sc.tl.rank_genes_groups(
                        adata_copy, groupby='_de_group',
                        groups=[cname], reference='Rest',
                        method=settings.DE_METHOD,
                    )
                    results[cname] = sc.get.rank_genes_groups_df(adata_copy, group=cname)
                except Exception as exc:
                    report(f"{cname} failed: {exc}")
        else:
            # Label each observation with the selected set it belongs to.
            # Sets can overlap (groups with each other, a group with a cluster):
            # observations in several of them are left out, so each set is
            # compared on what only it contains. Independent of what the map shows.
            labels = pd.Series(None, index=adata_copy.obs.index, dtype=object)
            n_member = np.zeros(len(labels), dtype=int)
            for cname in active:
                if cname in level_masks:
                    n_member += level_masks[cname].astype(bool)
                    labels[level_masks[cname].astype(bool)] = cname
            overlap = n_member > 1
            labels[overlap] = None
            n_overlap = int(overlap.sum())
            if n_overlap:
                report(f"Left out {n_overlap} {obs_word} in more than one selected {unit}")
            in_sel = labels.notna()
            sub = adata_copy[in_sel.values].copy()
            sub.obs['_de_group'] = labels[in_sel].astype(str).values
            tested_data = sub

            counts = sub.obs['_de_group'].value_counts()
            usable = [c for c in active if counts.get(c, 0) >= settings.DE_MIN_OBSERVATIONS]
            if len(usable) < 2:
                raise ValueError(f"\"Between selected\" needs at least 2 {units} "
                                 f"with ≥{settings.DE_MIN_OBSERVATIONS} {obs_word} each"
                                 + (" outside their overlap." if n_overlap else "."))

            report(f"Testing {len(usable)} {units} against each other")
            sc.tl.rank_genes_groups(
                sub, groupby='_de_group',
                groups=usable, reference='rest',
                method=settings.DE_METHOD,
            )
            for cname in usable:
                results[cname] = sc.get.rank_genes_groups_df(sub, group=cname)

        if not results:
            raise ValueError(f"No DE results for the selected {units} "
                             f"(each needs ≥{settings.DE_MIN_OBSERVATIONS} {obs_word}).")

        comparison_masks = {
            name: (np.asarray(tested_data.obs['_de_group'].astype(str) == name, dtype=bool)
                   if mode == 1 else np.asarray(level_masks[name], dtype=bool).copy())
            for name in results
        }
        # Snapshot observation metadata without duplicating the large expression matrix.
        expression = ad.AnnData(X=tested_data.X,
                                obs=pd.DataFrame(index=tested_data.obs_names.copy()),
                                var=tested_data.var.copy())
        return DEResult(results=results, expression=expression, masks=comparison_masks,
                        names=tuple(active), mode=mode, source=source, factor=factor,
                        n_overlap=n_overlap)

    def masks_at_level(self, masks, factor, level_obs_names):
        """Native-bin memberships -> memberships of the aggregated cells of `factor`.

        A cell belongs to a set when most of its native bins do (> 50%). Clusters
        computed at the same level map exactly: all bins of a cell share its label.
        """
        cell_names = self.sample._aggregated_cell_names[factor]   # native obs -> cell
        valid = cell_names.notna().to_numpy()
        codes = pd.Categorical(cell_names[valid], categories=level_obs_names).codes
        n_bins = np.bincount(codes, minlength=len(level_obs_names))
        level_masks = {}
        for name, mask in masks.items():
            in_set = np.bincount(codes, weights=np.asarray(mask, dtype=bool)[valid].astype(float),
                                 minlength=len(level_obs_names))
            level_masks[name] = in_set > 0.5 * np.maximum(n_bins, 1)
        report(f"Mapped {len(masks)} sets to {len(level_obs_names)} aggregated cells "
               f"(member when > 50% of their bins are)")
        return level_masks
