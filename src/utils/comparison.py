"""Coherence metrics between two samples of the same tissue at different resolutions.

The lower-resolution sample (s1, e.g. Visium spots) is the reference: each of its geometries is
compared with the (area-weighted) sum of the higher-resolution geometries (s2, e.g. Visium HD bins
or segmented cells) that fall inside it.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import shapely
from scipy import special
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

from .sample import Sample
from .progress import report


def compare_samples(
        s1,
        s2,
        min_coverage=0.8,
        predicate='contains',
        genes=None,
        n_hvg=2000,
        sample1_cluster_key=None,
        sample2_cluster_key=None,
        sample1_celltype_key=None,
        sample2_celltype_key=None,
        n_top_svg=100):
        
    """Compare a low-resolution and a high-resolution sample of the same tissue.

    Parameters
    ----------
    min_coverage : keep only low-res geometries whose area is covered more than this fraction
        by high-res geometries.
    predicate : spatial predicate used to match high-res geometries to low-res ones.
    genes : genes used for the HVG-restricted metrics. Defaults to `n_hvg` HVGs of the low-res sample.
    sample1_cluster_key, sample2_cluster_key : `adata.obs` columns with cluster labels.
    sample1_celltype_key : `s1.adata.obsm` key, or DataFrame indexed by obs names, with one column per
        cell type (deconvolution output, e.g. cell2location/RCTD).
    sample2_celltype_key : `s2.adata.obs` column with the cell type of each high-res geometry.

    Returns
    -------
    dict with
        "local": DataFrame, one row per matched low-res geometry, indexed like `s1.gdf`.
        "genes": DataFrame, one row per gene shared by both samples.
        "global": dict of scalar metrics.
    """
    if not isinstance(s1, Sample) or not isinstance(s2, Sample):
        raise ValueError("Both inputs must be instances of Sample.")

    if s1.gdf.geometry.area.mean() < s2.gdf.geometry.area.mean():
        s1, s2 = s2, s1

    # Geometric operations only involve the sequenced units (those present in the AnnData)
    gdf1, pos1 = _sequenced(s1)
    gdf2, pos2 = _sequenced(s2)

    report("Matching geometries...")
    pairs, coverage, overlap_area = _match_geometries(gdf1, gdf2, predicate, min_coverage)
    global_metrics = _tissue_metrics(gdf1, gdf2, overlap_area)
    del gdf2
    report(f"Number of low-res geometries after filtering by coverage: {pairs['left'].nunique()}")

    report("Building matched count matrices...")
    matched = _build_matched_matrices(s1, s2, gdf1, pos1, pos2, pairs) if not pairs.empty else None
    del pairs
    if matched is None:
        global_metrics["n_matched"] = 0
        return {"local": None, "genes": None, "global": global_metrics}

    L, H = matched["L"], matched["H"]
    report("Normalizing and selecting genes...")
    Ln, Hn = _log_normalize(L), _log_normalize(H)
    hvg = _select_genes(Ln, L, H, matched["genes"], genes, n_hvg)

    # ---- Local metrics ----
    report(f"Computing per-geometry metrics ({L.shape[0]} geometries)...")
    local = _profile_metrics(L, H, Ln, Hn, hvg)
    local.index = matched["labels"]
    local.insert(0, "obs_name", matched["obs_names"])
    local["coverage"] = coverage.loc[matched["spots"]].to_numpy()
    local["n_sample2"] = np.bincount(matched["rows"], minlength=len(matched["spots"]))
    local["distance_to_edge"] = _distance_to_edge(_tissue_hull(gdf1), matched["coords"])

    # ---- Per-gene metrics ----
    report(f"Computing per-gene metrics ({L.shape[1]} genes)...")
    gene_df = _gene_metrics(L, H, Ln, Hn, matched["genes"], hvg)

    report("Computing spatial autocorrelation...")
    W = _knn_weights(matched["coords"], k=6)
    gene_df.loc[gene_df.index[hvg], "morans_i_sample1"] = _morans_i(W, Ln[:, hvg].tocsc())
    gene_df.loc[gene_df.index[hvg], "morans_i_sample2"] = _morans_i(W, Hn[:, hvg].tocsc())

    # The count matrices are the largest objects: free them before anything else runs
    del L, H, Ln, Hn, matched["L"], matched["H"]

    # ---- Global metrics ----
    global_metrics.update({
        "n_matched": len(local),
        "n_sample1": len(gdf1),
        "n_sample2_matched": int(local["n_sample2"].sum()),
        "n_genes": len(gene_df),
        "n_hvg": len(hvg),
    })
    report("Computing global metrics...")
    global_metrics.update(_expression_metrics(local, gene_df))
    global_metrics.update(_spatial_metrics(gene_df, n_top_svg))
    global_metrics["local_metric_morans_i"] = _local_metric_autocorrelation(W, local)

    if sample1_cluster_key is not None and sample2_cluster_key is not None:
        report("Computing cluster agreement...")
        cluster_local, cluster_global = _cluster_metrics(s1, s2, matched, sample1_cluster_key, sample2_cluster_key)
        local = local.join(cluster_local)
        global_metrics.update(cluster_global)

    if sample1_celltype_key is not None and sample2_celltype_key is not None:
        report("Computing cell-type agreement...")
        celltype_local, celltype_global = _celltype_metrics(s1, s2, matched, sample1_celltype_key, sample2_celltype_key)
        local = local.join(celltype_local)
        global_metrics.update(celltype_global)

    return {"local": local, "genes": gene_df, "global": global_metrics}


# --------------------------------------------------------------------------------------------
# Matching
# --------------------------------------------------------------------------------------------

def _match_geometries(gdf1, gdf2, predicate, min_coverage):
    """(low-res row, high-res row) pairs, the fraction of each low-res geometry covered, and the
    total overlap area between both samples."""
    left = np.asarray(gdf1.geometry.values)
    right = np.asarray(gdf2.geometry.values)
    li, ri = shapely.STRtree(right).query(left, predicate='intersects')

    columns = ["left", "right", "intersection_area", "weight"]
    if len(li) == 0:
        return pd.DataFrame(columns=columns), pd.Series(dtype=float), 0.0

    # Only high-res geometries crossing a low-res boundary need an actual intersection
    shapely.prepare(left)
    inside = shapely.contains(left[li], right[ri])
    right_area = shapely.area(right[ri])
    intersection_area = right_area.copy()
    crossing = np.flatnonzero(~inside)
    intersection_area[crossing] = shapely.area(shapely.intersection(left[li[crossing]], right[ri[crossing]]))
    # Geometries of a sample do not overlap each other, so pairwise overlaps add up exactly
    overlap_area = float(intersection_area.sum())

    keep = inside if predicate == 'contains' else np.ones(len(li), dtype=bool)
    pairs = pd.DataFrame({
        "left": li[keep],
        "right": ri[keep],
        "intersection_area": intersection_area[keep],
        # Fraction of the high-res geometry inside the low-res one (1 with predicate='contains')
        "weight": intersection_area[keep] / right_area[keep],
    })

    # Keep only low-res geometries whose area is covered > min_coverage by high-res geometries
    covered = pairs.groupby("left")["intersection_area"].sum()
    coverage = covered / shapely.area(left[covered.index.to_numpy()])
    pairs = pairs[pairs["left"].isin(coverage.index[coverage > min_coverage])]
    return pairs, coverage, overlap_area


def _sequenced(sample):
    """Geometries of the sequenced units (rows of `sample.gdf` present in `sample.adata`) and the
    `adata` row of each of them."""
    obs = sample.gdf[sample.on].astype(str).to_numpy()
    pos = pd.Index(sample.adata.obs_names.astype(str)).get_indexer(obs)
    keep = np.flatnonzero(pos >= 0)
    return sample.gdf.iloc[keep], pos[keep]


def _build_matched_matrices(s1, s2, gdf1, pos1, pos2, pairs):
    """Low-res counts (L) and summed high-res counts (H) over the shared genes, one row per spot.
    `pairs` index rows of the sequenced geometries (`gdf1` and the one behind `pos2`)."""
    left = pairs["left"].to_numpy()
    right = pairs["right"].to_numpy()
    weight = pairs["weight"].to_numpy()

    spots = np.unique(left)
    rows = np.searchsorted(spots, left)
    A = sp.csr_matrix((weight.astype(np.float32), (rows, pos2[right])), shape=(len(spots), s2.adata.n_obs))

    shared_genes = s1.adata.var_names.intersection(s2.adata.var_names, sort=False)
    g1 = s1.adata.var_names.get_indexer(shared_genes)
    g2 = s2.adata.var_names.get_indexer(shared_genes)

    # float32 halves the memory of the (fairly dense) summed high-res matrix
    L = sp.csr_matrix(s1.adata.X[pos1[spots]])[:, g1].astype(np.float32)
    H = sp.csr_matrix(A @ s2.adata.X)[:, g2].astype(np.float32)

    # Spots with no counts on either side have undefined profiles
    nonempty = (np.asarray(L.sum(1)).ravel() > 0) & (np.asarray(H.sum(1)).ravel() > 0)
    if not nonempty.all():
        report(f"Dropping {(~nonempty).sum()} low-res geometries with zero counts in one of the samples")
        keep = np.flatnonzero(nonempty)
        remap = np.full(len(spots), -1)
        remap[keep] = np.arange(len(keep))
        in_kept = remap[rows] >= 0
        spots, L, H, A = spots[keep], L[keep], H[keep], A[keep]
        rows = remap[rows[in_kept]]
    if len(spots) == 0:
        return None

    # Genes detected in neither sample carry no information
    expressed = (L.getnnz(0) > 0) | (H.getnnz(0) > 0)
    L, H, shared_genes = L[:, expressed], H[:, expressed], shared_genes[expressed]

    centroids = gdf1.geometry.iloc[spots].centroid
    return {
        "spots": spots,
        "rows": rows,
        "labels": gdf1.index[spots],
        "obs_names": s1.adata.obs_names[pos1[spots]].astype(str),
        "coords": np.column_stack([centroids.x.to_numpy(), centroids.y.to_numpy()]),
        "A": A,
        "L": L,
        "H": H,
        "genes": shared_genes,
    }


# --------------------------------------------------------------------------------------------
# Normalisation and gene selection
# --------------------------------------------------------------------------------------------

def _log_normalize(M, target_sum=1e4):
    totals = np.asarray(M.sum(1, dtype=np.float64)).ravel()
    scale = np.divide(target_sum, totals, out=np.zeros_like(totals), where=totals > 0)
    return sp.csr_matrix(sp.diags(scale.astype(M.dtype)) @ M).log1p()


def _select_genes(Ln, L, H, shared_genes, genes, n_hvg):
    """Column indices of the genes used for HVG-restricted metrics."""
    if genes is not None:
        idx = shared_genes.get_indexer(pd.Index(genes))
        idx = idx[idx >= 0]
        if len(idx) == 0:
            raise ValueError("None of the requested genes is shared by both samples.")
        return np.sort(idx)

    # Only genes detected in both samples, otherwise they cannot agree by construction
    expressed = np.flatnonzero((L.getnnz(0) > 0) & (H.getnnz(0) > 0))
    n_top = min(n_hvg, len(expressed))
    X = Ln[:, expressed]
    try:
        import anndata as ad
        import scanpy as sc

        adata = ad.AnnData(X)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat')
        selected = np.flatnonzero(adata.var["highly_variable"].to_numpy())
    except Exception:
        _, var = _column_mean_var(X)
        selected = np.argsort(var)[::-1][:n_top]
    return np.sort(expressed[selected])


# --------------------------------------------------------------------------------------------
# Local (per low-res geometry) metrics
# --------------------------------------------------------------------------------------------

def _rowwise_pearson(a, b):
    a = a - a.mean(1, keepdims=True)
    b = b - b.mean(1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))


def _jensen_shannon(p, q):
    """Jensen-Shannon distance (base 2, in [0, 1]) between the rows of two composition matrices."""
    m = 0.5 * (p + q)
    divergence = 0.5 * (special.rel_entr(p, m).sum(1) + special.rel_entr(q, m).sum(1)) / np.log(2)
    return np.sqrt(np.clip(divergence, 0, None))


def _profile_metrics(L, H, Ln, Hn, hvg, chunk_size=256):
    total_l = np.asarray(L.sum(1, dtype=np.float64)).ravel()
    total_h = np.asarray(H.sum(1, dtype=np.float64)).ravel()
    detected_l, detected_h = L.getnnz(1), H.getnnz(1)

    Ln_hvg, Hn_hvg = Ln[:, hvg], Hn[:, hvg]
    pearson, js = [], []
    for start in range(0, L.shape[0], chunk_size):
        sl = slice(start, start + chunk_size)
        pearson.append(_rowwise_pearson(Ln_hvg[sl].toarray(), Hn_hvg[sl].toarray()))
        l, h = L[sl].toarray(), H[sl].toarray()
        js.append(_jensen_shannon(l / total_l[sl, None], h / total_h[sl, None]))

    return pd.DataFrame({
        "total_counts_sample1": total_l,
        "total_counts_sample2": total_h,
        "log2_count_ratio": np.log2(total_h / total_l),
        "detected_genes_sample1": detected_l,
        "detected_genes_sample2": detected_h,
        # Expression profile agreement on the informative genes
        "pearson_hvg": np.concatenate(pearson),
        # Composition divergence over all genes
        "jensen_shannon": np.concatenate(js),
    })


def _label_composition(A, labels):
    """Weighted count of every high-res label inside each low-res geometry."""
    codes, categories = pd.factorize(pd.Series(labels))
    valid = np.flatnonzero(codes >= 0)
    onehot = sp.csr_matrix((np.ones(len(valid)), (valid, codes[valid])),
                           shape=(len(labels), len(categories)))
    return pd.DataFrame((A @ onehot).toarray(), columns=pd.Index(categories).astype(str))


def _cluster_metrics(s1, s2, matched, sample1_key, sample2_key):
    composition = _label_composition(matched["A"], s2.adata.obs[sample2_key].to_numpy())
    totals = composition.sum(1)
    has_labels = totals > 0

    majority = composition.idxmax(1).where(has_labels)
    purity = (composition.max(1) / totals).where(has_labels)
    sample1_labels = s1.adata.obs[sample1_key].set_axis(s1.adata.obs_names.astype(str))
    sample1 = pd.Series(sample1_labels.reindex(matched["obs_names"]).to_numpy(), dtype=object)

    local = pd.DataFrame({
        "cluster_sample1": sample1.to_numpy(),
        "cluster_sample2_majority": majority.to_numpy(),
        "cluster_sample2_purity": purity.to_numpy(),
        # Only meaningful when both samples share the label space (e.g. transferred labels)
        "cluster_match": (sample1.astype(str) == majority.astype(str)).to_numpy(),
    }, index=matched["labels"])

    both = sample1.notna().to_numpy() & majority.notna().to_numpy()
    a, b = sample1[both].astype(str), majority[both].astype(str)
    return local, {
        "cluster_ari": adjusted_rand_score(a, b) if both.any() else np.nan,
        "cluster_match_rate": float((a.to_numpy() == b.to_numpy()).mean()) if both.any() else np.nan,
        "cluster_sample2_purity_mean": float(purity.mean()),
    }


def _celltype_metrics(s1, s2, matched, sample1_celltype_key, sample2_key):
    sample1 = s1.adata.obsm[sample1_celltype_key] if isinstance(sample1_celltype_key, str) else sample1_celltype_key
    if not isinstance(sample1, pd.DataFrame):
        raise ValueError("sample1_celltype_key must be a DataFrame (or obsm key holding one) with one column per cell type.")
    sample1 = sample1.copy()
    sample1.index = sample1.index.astype(str)
    sample1.columns = sample1.columns.astype(str)

    sample2 = _label_composition(matched["A"], s2.adata.obs[sample2_key].to_numpy())
    celltypes = sample1.columns.intersection(sample2.columns, sort=False)
    if len(celltypes) == 0:
        raise ValueError("No cell type names shared between sample1_celltype_key and sample2_celltype_key.")

    p = sample1.reindex(matched["obs_names"])[celltypes].fillna(0).to_numpy(dtype=float)
    q = sample2[celltypes].to_numpy(dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        p = p / p.sum(1, keepdims=True)
        q = q / q.sum(1, keepdims=True)

    local = pd.DataFrame({
        "celltype_jensen_shannon": _jensen_shannon(p, q),
    }, index=matched["labels"])

    per_celltype = {ct: _pearson(p[:, i], q[:, i]) for i, ct in enumerate(celltypes)}
    return local, {
        "celltype_jensen_shannon_mean": float(np.nanmean(local["celltype_jensen_shannon"])),
        "celltype_spatial_pearson": per_celltype,
    }


# --------------------------------------------------------------------------------------------
# Per-gene and global expression metrics
# --------------------------------------------------------------------------------------------

def _pearson(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def _column_mean_var(X):
    mean = np.asarray(X.mean(0, dtype=np.float64)).ravel()
    var = np.asarray(X.multiply(X).mean(0, dtype=np.float64)).ravel() - mean ** 2
    return mean, np.clip(var, 0, None)


def _columnwise_pearson(X, Y):
    mx, vx = _column_mean_var(X)
    my, vy = _column_mean_var(Y)
    cov = np.asarray(X.multiply(Y).mean(0, dtype=np.float64)).ravel() - mx * my
    with np.errstate(invalid='ignore', divide='ignore'):
        return cov / np.sqrt(vx * vy)


def _gene_metrics(L, H, Ln, Hn, genes, hvg):
    n = L.shape[0]
    total_l = np.asarray(L.sum(0, dtype=np.float64)).ravel()
    total_h = np.asarray(H.sum(0, dtype=np.float64)).ravel()

    return pd.DataFrame({
        "mean_sample1": total_l / n,
        "mean_sample2": total_h / n,
        "capture_ratio": np.divide(total_h, total_l, out=np.full(len(genes), np.nan), where=total_l > 0),
        "zero_fraction_sample1": 1 - L.getnnz(0) / n,
        "zero_fraction_sample2": 1 - H.getnnz(0) / n,
        # Does the gene have the same spatial pattern across matched spots?
        "spatial_pearson": _columnwise_pearson(Ln, Hn),
        "is_hvg": np.isin(np.arange(len(genes)), hvg),
        "morans_i_sample1": np.nan,
        "morans_i_sample2": np.nan,
    }, index=pd.Index(genes, name="gene"))


def _expression_metrics(local, gene_df):
    total_l = local["total_counts_sample1"].to_numpy()
    total_h = local["total_counts_sample2"].to_numpy()

    pb_l = gene_df["mean_sample1"].to_numpy()
    pb_h = gene_df["mean_sample2"].to_numpy()
    pb_ln = np.log1p(pb_l / pb_l.sum() * 1e6)
    pb_hn = np.log1p(pb_h / pb_h.sum() * 1e6)

    # Genes detected in at least 10% of matched spots in both samples
    detected = (gene_df["zero_fraction_sample1"] < 0.9) & (gene_df["zero_fraction_sample2"] < 0.9)

    return {
        "pseudobulk_pearson": _pearson(pb_ln, pb_hn),
        "capture_ratio": float(total_h.sum() / total_l.sum()),
        "total_counts_pearson": _pearson(np.log1p(total_l), np.log1p(total_h)),
        "pearson_hvg_median": float(np.nanmedian(local["pearson_hvg"])),
        "jensen_shannon_median": float(np.nanmedian(local["jensen_shannon"])),
        "gene_spatial_pearson_median": float(gene_df.loc[detected, "spatial_pearson"].median()),
        "gene_spatial_pearson_median_hvg": float(gene_df.loc[gene_df["is_hvg"], "spatial_pearson"].median()),
    }


# --------------------------------------------------------------------------------------------
# Spatial autocorrelation
# --------------------------------------------------------------------------------------------

def _knn_weights(coords, k=6):
    """Row-standardised k-nearest-neighbour spatial weights."""
    n = len(coords)
    k = min(k, n - 1)
    if k < 1:
        return None
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(coords)
    rows = np.repeat(np.arange(n), k)
    return sp.csr_matrix((np.full(n * k, 1.0 / k), (rows, idx[:, 1:].ravel())), shape=(n, n))


def _morans_i(W, X, chunk_size=256):
    """Moran's I of every column of X (row-standardised W, so S0 = n), densified in column chunks."""
    if W is None:
        return np.full(X.shape[1], np.nan)
    out = []
    for start in range(0, X.shape[1], chunk_size):
        Z = X[:, start:start + chunk_size]
        Z = Z.toarray() if sp.issparse(Z) else np.asarray(Z, dtype=float)
        Z = Z - Z.mean(0)
        with np.errstate(invalid='ignore', divide='ignore'):
            out.append((Z * (W @ Z)).sum(0) / (Z ** 2).sum(0))
    return np.concatenate(out)


def _spatial_metrics(gene_df, n_top_svg):
    svg = gene_df.dropna(subset=["morans_i_sample1", "morans_i_sample2"])
    k = min(n_top_svg, len(svg))
    top_l = set(svg["morans_i_sample1"].nlargest(k).index)
    top_h = set(svg["morans_i_sample2"].nlargest(k).index)
    return {
        "morans_i_pearson": _pearson(svg["morans_i_sample1"].to_numpy(), svg["morans_i_sample2"].to_numpy()),
        "svg_top_jaccard": len(top_l & top_h) / len(top_l | top_h) if k else np.nan,
    }


def _local_metric_autocorrelation(W, local):
    """Moran's I of the main local metric maps: high values mean the disagreement is regional."""
    numeric = local[["pearson_hvg", "jensen_shannon", "log2_count_ratio"]].replace([np.inf, -np.inf], np.nan)
    numeric = numeric.fillna(numeric.mean())
    values = _morans_i(W, numeric.to_numpy())
    return dict(zip(numeric.columns, values.astype(float).tolist()))


# --------------------------------------------------------------------------------------------
# Tissue geometry
# --------------------------------------------------------------------------------------------

def _tissue_hull(gdf, ratio=0.1):
    """Tissue shape used only for distances to the tissue edge: concave hull of the geometry
    centroids (never enlarges the tissue beyond its outermost geometries)."""
    centroids = gdf.geometry.centroid
    return shapely.concave_hull(shapely.multipoints(shapely.points(
        np.column_stack([centroids.x.to_numpy(), centroids.y.to_numpy()]))), ratio=ratio)


def _tissue_metrics(gdf1, gdf2, overlap_area):
    """Tissue overlap of the two clouds of geometries (no buffering: gaps are not tissue)."""
    area1 = float(gdf1.geometry.area.sum())
    area2 = float(gdf2.geometry.area.sum())
    return {
        "tissue_area_sample1": area1,
        "tissue_area_sample2": area2,
        "tissue_intersection_area": overlap_area,
        "tissue_dice": 2 * overlap_area / (area1 + area2) if area1 + area2 else np.nan,
    }


def _distance_to_edge(outline, coords):
    """Distance from each point to the tissue boundary (0 outside the tissue)."""
    points = shapely.points(coords)
    distance = shapely.distance(points, outline.boundary)
    return np.where(shapely.contains(outline, points), distance, 0.0)
