"""In-app help: the startup tutorial and the "?" buttons across the platform.

All user-facing explanations live here so they stay consistent and easy to edit.
Each entry describes what the code actually does; update it with the feature.
"""

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
        It is loaded once and reused in every layer's Deconv. (deconvolution) tab.<br>
        Examples use their configured local reference, if available; otherwise they run without one.
        Reference uploads are available only with your own data.<br>
        The <b>Human CRC</b> example loads its Chromium counts, single-cell CSV annotations and spatial
        Parquet metadata directly on the server, with no upload. Deconvolution defaults to the <b>Level1</b> cell-type labels;
        you can select Level2 or another reference annotation column."""),
    "preprocessing": ("What to compute", """
        <b>Bins</b> are always loaded.<br>
        <b>Pseudospots</b> group the HD bins into Visium v2-like 55 µm spots, e.g. to compare with
        standard Visium data.<br>
        <b>Cell segmentation</b>: each method you tick finds cells (or nuclei) in the tissue image and
        assigns the bins to them; every method becomes its own cells layer, so they can be compared.
        Results are cached, so loading the same image with the same settings is instant.
        Methods run in their own containers; their parameters are shown when ticked.
        Available methods depend on the server's setup. Each method is a plugin: developers can add
        their own (<a href='static/assets/tutorial.html#app-plugins' target='_blank' rel='noopener'>how</a>)."""),
    "progress": ("Loading progress", """
        One bar per task. Percentages count completed steps, not time remaining. A failed
        segmentation method turns its bar red with the reason; the rest of the data still loads.
        If a problem persists, report it with the message shown at
        <a href='https://github.com/LuisAlonsoEsteban/VisiumVisualizer/issues' target='_blank' rel='noopener'>GitHub issues</a>
        or to <a href='mailto:laesteban@unav.es'>laesteban@unav.es</a>."""),

    # ------------------------------------------------------------- layer tabs
    "general": ("Showing a layer", """
        Show or hide the whole layer with the eye in the layer list. Only observations with counts are
        loaded: bins, spots and cells without counts are left out.<br>
        <b>Opacity</b>: maximum opacity. When coloring by a gene, opacity follows its expression up to
        this value; with categories (clusters, groups) every observation uses it.<br>
        <b>Clear selection</b> deselects the layer's observations; <b>Download data</b> saves the layer as an
        AnnData file (.h5ad) with its clusters, groups and deconvolution results."""),
    "color": ("Coloring", """
        Choose a <b>gene or metadata column</b> to display:
        <ul>
          <li>Genes and numeric columns: choose one color with the color picker; opacity follows the value.</li>
          <li><b>Expression values</b>: show genes as raw counts, normalized counts (per 10,000 counts of
              each observation) or log-normalized counts (log1p). Only the display changes; analyses
              and group conditions use counts.</li>
          <li>Categorical columns (clusters, groups, ...): the single-color picker is disabled. Each category has its own color, the same colors used
              in the cluster/group lists and the DE results.</li>
        </ul>"""),
    "clustering": ("Clustering", """
        Choose <b>Leiden, K-means, DBSCAN, HDBSCAN, Gaussian mixtures or agglomerative clustering</b>.
        All methods use normalized, log-transformed highly variable gene expression and PCA.
        Leiden uses a neighbor graph and resolution; K-means, mixtures and agglomerative clustering
        use a chosen number of clusters. DBSCAN uses an epsilon distance in PCA space and minimum
        samples; HDBSCAN uses minimum cluster size and minimum samples. Density-method outliers
        appear as <b>Noise</b>. Mixtures offer covariance types and agglomerative clustering offers linkage choices.
        Clusters appear first; marker genes are ranked afterwards. Noise and groups with fewer than
        two observations have no marker ranking. A single resulting group also has no marker ranking.
        For bins, clustering runs at the <b>displayed aggregation level</b>. While its clusters are shown,
        that level is the minimum: dynamic zoom and coarser fixed levels remain available.
        Coarser bins show their most common cluster label. Tick or untick clusters to show or hide them."""),
    "deconvolution": ("Deconvolution", """
        The single-cell <b>.h5ad or .h5 reference</b> uploaded on the loading screen is reused
        for every layer. Its .obs columns (including an optional annotation CSV or Parquet file) appear here. Select its
        cell-type column and run deconvolution. Both datasets must use matching, unique gene identifiers and linear,
        nonnegative expression in <code>.X</code> (not log-transformed values).<br>
        Choose <b>NNLS, spacexr (RCTD), cell2location, MuSiC, or CIBERSORTx</b> and set its parameters;
        methods not installed on this server are not listed. Each method is a plugin, and developers can
        add their own (<a href='static/assets/tutorial.html#app-plugins' target='_blank' rel='noopener'>how</a>).
        NNLS estimates expression mixture fractions rather than calibrated cell counts.
        spacexr (RCTD) needs raw integer counts; set <b>RCTD mode</b> (doublet by default, full or multi)
        and <b>Minimum target UMIs</b> (default 100), plus normalization cutoffs and reference filters.
        It can take hours on a whole 8 µm sample. <b>Cancel</b> stops a running deconvolution; no result is saved.
        cell2location needs an expected number of cells per location and offers posterior mean abundance
        as well as relative fractions. MuSiC requires a donor/sample column with multiple subjects.
        CIBERSORTx generates a signature from the reference; it is offered only where the server operator has
        installed it with their own CIBERSORTx credentials.<br>
        Runs continue in the background only while this browser session stays open.<br>
        Results cover all native observations of this layer, independently of the viewport.
        Empty observations stay unassigned. Select a saved run, readout, and cell type to show it on the map.
        Coarse bin levels average fractions and sum expected cell abundance. Previous runs and their
        parameters are stored in the layer's AnnData."""),
    "rois": ("ROIs", """
        Save named regions in this layer for the current session. Choose <b>Bounding box</b>,
        <b>Lasso</b>, or <b>Seed + radius</b>, then click <b>Draw</b> and draw on the map.
        A seed radius uses µm when calibrated, otherwise full-resolution pixels. Radii and ruler
        measurements account for image scale and padding.
        On bins layers, <b>Find barcode / bin index</b> searches native barcodes with autocomplete.
        Select an exact barcode to preview and focus its radius ROI. The seed is the native bin's
        center, even outside the viewport or at a coarser display level. Changing the radius updates
        the draft; click <b>Preview barcode</b> to rebuild and focus it, then <b>Save ROI</b>.
        On bins layers, choose <b>All bins in a group</b> as the seed source and select a group.
        Every native bin in that group supplies a seed; all radius neighborhoods are combined
        into <b>one ROI</b>. Choose <b>Full ROI</b> to keep the complete neighborhood, or
        <b>ROI minus group</b> to subtract the source-group bin shapes and exclude its members.
        The preview updates when the group, radius, or ROI region changes, including
        disconnected regions and holes. Complex outlines use an overview image to keep the map
        responsive; zoomed edges can look pixelated, while membership uses the exact region.
        Inspect the draft and click <b>Save ROI</b> to keep it; Cancel discards it.<br>
        Counts include all sequenced native observations whose shapes intersect the ROI, including
        boundary contact, regardless of the viewport or displayed bin resolution.<br>
        <b>Focus</b> fits the viewport, <b>Select</b> links the selection to the map and embeddings,
        and <b>Create group</b> creates a snapshot for DE. Empty regions can be focused but cannot
        create groups. <b>Delete</b> also deletes the ROI's linked group and clears affected DE results.
        Outlines can be hidden independently. ROI shapes are fixed after saving; redraw to replace them."""),
    "groups": ("Groups", """
        Define your own groups from conditions on genes or metadata, e.g.
        <code>EPCAM &gt; 2</code> <b>AND</b> <code>COL1A1 &gt; 0</code> (one condition at a time, combined with
        AND/OR). Categorical columns accept <code>==</code> and a category name.
        Build the conditions, name the group and create it.<br>
        An observation in several groups is drawn with the group created last. Tick or untick groups to
        show or hide them; this does not change differential expression.<br>
        <b>Combine groups</b> makes a new group from two existing ones: their <b>intersection</b> (A ∩ B),
        <b>union</b> (A ∪ B) or <b>difference</b> (A − B), or the <b>complement</b> of one (not A: every
        sequenced observation outside it). <i>Take the opposite set</i> keeps everything outside the result.
        The size is shown before creating it, and its members stay fixed if A or B change later."""),
    "de": ("Differential expression", """
        Finds the genes that differ between groups and/or clusters (Wilcoxon test on log1p library-size-normalized counts).
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
    "metadata_comparison": ("Metadata comparison", """
        Compare <b>numeric metadata</b> between two or more groups and/or clusters.
        Select <b>The rest</b> to compare one or more selected sets against observations outside
        <b>all defined custom groups</b> (including hidden and unselected groups) and outside the
        compared set itself. Each selected group/cluster gets its own plot against its own rest.
        If nothing is left outside the groups and the set, The rest is empty and cannot be compared.
        Choose a saved deconvolution run to compare its cell-type fractions or abundance,
        or select arbitrary numeric .obs columns. The default summary is the <b>mean</b>;
        median is also available.<br>
        Results appear <b>below the map</b>: horizontal bars put cell types or columns on the y-axis,
        values on the x-axis, and give each group its own bar and color. Hover shows the value
        and number of finite observations.
        Comparisons use native observations regardless of viewport or displayed bin aggregation.
        Missing and infinite values are omitted per column. Overlapping observations contribute
        to each group. Fractions are averaged, without renormalizing the selected cell types.
        These are descriptive summaries, without significance tests. DE remains a separate tab."""),
    "parameters": ("Parameters", """
        Extra drawing parameters of this layer: line color, line width, fill and line opacity.
        Each is either a <b>fixed value</b> or taken <b>from data</b> (a gene or numeric column scaled to
        an output range; for colors, a categorical column). <b>Apply</b> redraws the map at the same view;
        <b>Reset</b> returns the attribute to its default. Color and opacity themselves are set in the
        Display tab."""),

    # ----------------------------------------------------------- other panels
    "bin_resolution": ("Bin resolution", """
        <b>Dynamic</b>: zoomed out, bins are summed into larger cells (2 × 2 up to 12 × 12) to keep the
        map fast; zooming in shows native bins again.<br>
        <b>Fixed</b>: always show one level. While bin clusters are shown, their computed level is the
        minimum. Dynamic mode follows zoom above that minimum; finer levels are hidden from the selector."""),
    "embedding": ("Embeddings", """
        PCA, UMAP or t-SNE of a layer, computed on the same preprocessing as clustering (so a UMAP after
        clustering only computes the UMAP itself). For bins it uses the displayed aggregation level.
        Select points to highlight the same observations on the map, and the other way round.
        <b>Color by</b> colors the points by a gene (on the layer's expression scale) or a metadata column;
        clusters and groups keep their assigned colors, and a new clustering is shown automatically."""),
    "metadata": ("Summary", """
        Per layer: quality metrics (counts, detected genes and % mitochondrial reads per observation), the marker genes of the last
        clustering, and differential-expression results."""),
}


# ------------------------------------------------------------- tutorial page
# Standalone tutorial pages in hd_viewer/static, served at static/assets/ (relative,
# so they work under any URL prefix)
TUTORIAL_URL = "static/assets/tutorial.html"
BIOLOGICAL_USE_CASE_URL = "static/assets/biological-use-case.html"

# Section of the tutorial each help topic links to
TUTORIAL_ANCHORS = {
    "data_source": "load-files", "preprocessing": "load-preprocessing", "progress": "load-progress",
    "general": "app-general", "color": "app-color", "clustering": "app-clustering",
    "groups": "app-groups", "de": "app-de", "parameters": "app-parameters",
    "rois": "app-rois", "metadata_comparison": "app-metadata-comparison",
    "deconvolution": "app-deconvolution",
    "bin_resolution": "app-bin-resolution", "embedding": "app-embedding", "metadata": "app-summary",
}


def help_topics():
    """Every help topic as {topic: {title, html, tutorial}} (tutorial: a URL or None)."""
    return {topic: dict(title=title, html=" ".join(body.split()),
                        tutorial=f"{TUTORIAL_URL}#{TUTORIAL_ANCHORS[topic]}" if topic in TUTORIAL_ANCHORS else None)
            for topic, (title, body) in HELP.items()}
