"""In-app help: the startup tutorial and the "?" buttons across the platform.

All user-facing explanations live here so they stay consistent and easy to edit.
Each entry describes what the code actually does; update it with the feature.
"""

from bokeh.layouts import row
from bokeh.models import Div, HelpButton, Tooltip
from bokeh.models.dom import HTML

# topic -> (title, HTML body)
HELP = {
    # ---------------------------------------------------------------- startup
    "data_source": ("Your data", """
        Upload the Space Ranger output of one Visium HD sample:
        <ol>
          <li><b>Expression (.h5)</b>: <code>binned_outputs/square_XXXum/filtered_feature_bc_matrix.h5</code></li>
          <li><b>Positions (.parquet)</b>: <code>binned_outputs/square_XXXum/spatial/tissue_positions.parquet</code>,
              from the same bin size</li>
          <li><b>Resolution</b>: that bin size (2, 8 or 16 µm)</li>
          <li><b>Tissue image</b>: the full-resolution image given to Space Ranger (.btf, .tif or .png)</li>
        </ol>
        Or pick a ready-made <b>example</b> dataset.<br>
        When uploading your own spatial data, optionally upload a single-cell <b>.h5ad or .h5 reference</b>.
        AnnData .obs annotations load automatically; for 10x count matrices, upload a CSV or Parquet file with
        cell barcodes in its first column and cell-type annotations in the other columns.
        Parquet can also store barcodes in its index or a barcode column.
        It is loaded once and reused in every layer's Deconvolution tab.<br>
        Examples use their configured local reference, if available; otherwise they run without one.
        Reference uploads are available only with your own data.<br>
        The <b>Human CRC</b> example loads its Chromium counts, single-cell CSV annotations and spatial
        Parquet metadata directly from disk. Deconvolution defaults to the <b>Level1</b> cell-type labels;
        you can select Level2 or another reference annotation column."""),
    "preprocessing": ("What to compute", """
        <b>Bins</b> are always loaded.<br>
        <b>Pseudospots</b> group the HD bins into Visium v2-like 55 µm spots, e.g. to compare with
        standard Visium data.<br>
        <b>Cell segmentation</b>: each method you tick finds cells (or nuclei) in the tissue image and
        assigns the bins to them; every method becomes its own cells layer, so they can be compared.
        Results are cached, so loading the same image with the same settings is instant.
        Methods run in their own containers; their parameters are shown when ticked."""),
    "progress": ("Loading progress", """
        One bar per task. Percentages count completed steps, not time remaining. A failed
        segmentation method turns its bar red with the reason; the rest of the data still loads.
        The full error is printed in the server console and kept in the method's
        <code>run.log</code>."""),

    # ------------------------------------------------------------- layer tabs
    "general": ("General", """
        <b>Show</b> or hide the whole layer, and the <b>not-sequenced</b> observations (those without
        counts).<br>
        <b>Alpha</b>: maximum opacity. When colouring by a gene, opacity follows its expression up to
        this value; with categories (clusters, groups) every observation uses it."""),
    "color": ("Color", """
        Choose a <b>gene or metadata column</b> to display:
        <ul>
          <li>Genes and numeric columns: one colour (the colour picker) whose opacity follows the value.</li>
          <li>Categorical columns (clusters, groups, ...): one colour per category, the same colours used
              in the cluster/group lists and the DE results.</li>
        </ul>"""),
    "clustering": ("Clustering", """
        Choose <b>Leiden, K-means, DBSCAN, HDBSCAN, Gaussian mixtures or agglomerative clustering</b>.
        All methods use normalised, log-transformed highly variable gene expression and PCA.
        Leiden uses a neighbour graph and resolution; K-means, mixtures and agglomerative clustering
        use a chosen number of clusters. DBSCAN uses an epsilon distance in PCA space and minimum
        samples; HDBSCAN uses minimum cluster size and minimum samples. Density-method outliers
        appear as <b>Noise</b>. Mixtures offer covariance types and agglomerative clustering offers linkage choices.
        Clusters appear first; marker genes are ranked afterwards. Noise and groups with fewer than
        two observations have no marker ranking. A single resulting group also has no marker ranking.
        For bins, clustering runs at the <b>displayed aggregation level</b> and keeps the view at that
        level while the clusters are shown. Tick or untick clusters to show or hide them."""),
    "deconvolution": ("Deconvolution", """
        The single-cell <b>.h5ad or .h5 reference</b> uploaded on the loading screen is reused
        for every layer. Its .obs columns (including an optional annotation CSV or Parquet file) appear here. Select its
        cell-type column and run deconvolution. Both datasets must use matching, unique gene identifiers and linear,
        nonnegative expression in <code>.X</code> (not log-transformed values).<br>
        Choose <b>NNLS, cell2location, MuSiC, or CIBERSORTx</b> and set its parameters.
        NNLS estimates expression mixture fractions rather than calibrated cell counts.
        cell2location needs an expected number of cells per location and offers posterior mean abundance
        as well as relative fractions. MuSiC requires a donor/sample column with multiple subjects.
        CIBERSORTx generates a signature from the reference and requires a configured container and credentials.<br>
        Results cover all native observations of this layer, independently of the viewport.
        Empty observations stay unassigned. Select a saved run, readout, and cell type to show it on the map.
        Coarse bin levels average fractions and sum expected cell abundance. Previous runs and their
        parameters are stored in the layer's AnnData."""),
    "rois": ("ROIs", """
        Save named regions in this layer for the current session. Choose <b>Bounding box</b>,
        <b>Lasso</b>, or <b>Seed + radius</b>, then click <b>Draw</b> and draw on the map.
        A seed radius uses µm when calibrated, otherwise full-resolution pixels. Radii and ruler
        measurements account for image scale and padding.
        On bins layers, choose <b>All bins in a group</b> as the seed source and select a group.
        Every native bin in that group supplies a seed; all radius neighbourhoods are combined
        into <b>one ROI</b>. Choose <b>Full ROI</b> to keep the complete neighbourhood, or
        <b>ROI minus group</b> to subtract the source-group bin shapes and exclude its members.
        The preview updates when the group, radius, or ROI region changes, including
        disconnected regions and holes. Complex outlines use an overview image to keep the map
        responsive; zoomed edges can look pixelated, while membership uses the exact region.
        Click Save ROI to keep it.
        Inspect the draft and click <b>Save ROI</b>; Cancel discards it.<br>
        Counts include all sequenced native observations whose shapes intersect the ROI, including
        boundary contact, regardless of the viewport or displayed bin resolution.<br>
        <b>Focus</b> fits the viewport, <b>Select</b> links the selection to the map and embeddings,
        and <b>Create group</b> creates a snapshot for DE. Empty regions can be focused but cannot
        create groups. <b>Delete</b> also deletes the ROI's linked group and clears affected DE results.
        Outlines can be hidden independently. ROI shapes are fixed after saving; redraw to replace them."""),
    "groups": ("Groups", """
        Define your own groups from conditions on genes or metadata, e.g.
        <code>Ttr &gt; 2 AND Pmch &gt; 0</code>. Build the conditions, name the group and create it.<br>
        An observation in several groups is drawn with the group created last. Tick or untick groups to
        show or hide them; this does not change differential expression."""),
    "de": ("Differential expression", """
        Finds the genes that differ between groups and/or clusters (Wilcoxon test on log1p counts).
        <ul>
          <li><b>vs Rest</b>: each ticked group or cluster against everything else.</li>
          <li><b>Between selected</b>: the ticked ones against each other, e.g. a group vs a cluster.
              Observations in more than one of them are left out.</li>
        </ul>
        For bins, <b>Resolution</b> tests native bins or aggregated cells (summed counts); a cell belongs
        to a group or cluster when more than half of its bins do.<br>
        Results appear below the map: click a column header to sort, and search any gene.
        Plots support hover, pan, zoom and selection. Hover dots or violins for gene statistics;
        click a legend entry to hide or show a comparison."""),
    "parameters": ("Parameters", """
        Extra drawing parameters of this layer: line colour, line width, fill and line opacity.
        Each is either a <b>fixed value</b> or taken <b>from data</b> (a gene or numeric column scaled to
        an output range; for colours, a categorical column). Edit (✎) or remove (✕) them from the list.
        Colour and alpha themselves are set in the Color and General tabs."""),

    # ----------------------------------------------------------- other panels
    "bin_resolution": ("Bin resolution", """
        <b>Dynamic</b>: zoomed out, bins are summed into larger cells (2 × 2 up to 12 × 12) to keep the
        map fast; zooming in shows native bins again.<br>
        <b>Fixed</b>: always show one level. While bin clusters are shown, the level they were computed
        at is used."""),
    "embedding": ("Embeddings", """
        PCA, UMAP or t-SNE of a layer, computed on the same preprocessing as clustering (so a UMAP after
        clustering only computes the UMAP itself). For bins it uses the displayed aggregation level.
        Select points to highlight the same observations on the map, and the other way round."""),
    "metadata": ("Summary", """
        Per layer: quality metrics (counts, detected genes and % mitochondrial reads per observation), the marker genes of the last
        clustering, and differential-expression results."""),
}


def help_html(topic):
    title, body = HELP.get(topic, ("",""))
    anchor = TUTORIAL_ANCHORS.get(topic, "")
    more = (f"<div style='margin-top:6px'><a href='{TUTORIAL_URL}#{anchor}' target='_blank' "
            f"rel='noopener'>Read more in the tutorial →</a></div>") if anchor else ""
    return (f"<div style='max-width:340px;font-size:12px;line-height:1.45'>"
            f"<b style='font-size:13px'>{title}</b><div style='margin-top:4px'>{body}</div>{more}</div>")


def help_button(topic, position="right"):
    """Small "?" button explaining `topic` on hover or click."""
    return HelpButton(tooltip=Tooltip(content=HTML(help_html(topic)), position=position),
                      align="center", margin=(0, 0, 0, 6))


def tutorial_link():
    """Link opening the tutorial in a new browser tab."""
    return Div(text=f"<a href='{TUTORIAL_URL}' target='_blank' rel='noopener' "
                    f"style='color:#88f5dc;font-size:12px;white-space:nowrap'>📖 Tutorial</a>",
               align="center", margin=(0, 0, 0, 12))


def help_header(text, topic, style="color:#888888"):
    """A section title with its help button."""
    return row(Div(text=f"<b style='{style}'>{text}</b>"), help_button(topic), margin=0)


# ------------------------------------------------------------- tutorial page
# Standalone tutorial page in bokeh_server/static (served next to the app)
TUTORIAL_URL = "/hd_viewer/static/tutorial.html"

# Section of the tutorial each help topic links to
TUTORIAL_ANCHORS = {
    "data_source": "load-files", "preprocessing": "load-preprocessing", "progress": "load-progress",
    "general": "app-general", "color": "app-color", "clustering": "app-clustering",
    "groups": "app-groups", "de": "app-de", "parameters": "app-parameters",
    "rois": "app-rois",
    "deconvolution": "app-deconvolution",
    "bin_resolution": "app-bin-resolution", "embedding": "app-embedding", "metadata": "app-summary",
}

TUTORIAL_CARD_HTML = f"""
<div style="display:flex;flex-wrap:wrap;align-items:center;gap:12px 20px">
  <div style="flex:1 1 320px">
    <h2 style="margin:0 0 4px;font-size:19px">New here?</h2>
    <p style="margin:0;color:#b7cbd1">The tutorial explains every part of the platform: which files to load,
       what each option computes, and how to explore the results.</p>
  </div>
  <a href="{TUTORIAL_URL}" target="_blank" rel="noopener"
     style="background:#5ee0c2;color:#102a24;font-weight:600;text-decoration:none;
            padding:12px 18px;border-radius:10px;white-space:nowrap">📖 Open the tutorial</a>
</div>
"""
