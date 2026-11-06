# React interface: parity checklist (milestones 1–6)

Compared with `bokeh serve hd_viewer` on the same dataset. "Verified" means
covered by an automated test (tests/test_api.py, tests/test_core_*.py,
frontend Vitest/Playwright) or checked in a browser.

| Feature | Bokeh UI | React UI | Status |
|---|---|---|---|
| Landing page, tutorial and walkthrough links, citation, contact | yes | yes | verified (Playwright) |
| Load an example dataset, per-task progress | yes | yes (events socket, polling fallback) | verified (API + Playwright, 16 µm mouse brain) |
| Upload files in chunks (validation, size limit, replace/cancel) | yes | yes (HTTP PUT) | verified (API tests); browser upload not yet automated |
| Bin size of uploads, optional single-cell reference | yes | yes | API-tested |
| Pseudospots option | yes | yes | API path shared with core pipeline |
| Segmentation methods and parameter forms from manifests | yes | yes | forms render (Playwright); container runs not testable here |
| Tissue image tiles (WebP, immutable caching) | yes | yes (`/api/tiles`, one port) | verified |
| Lazy loading of bins with aggregation levels and hysteresis | yes | yes (shared core.viewport / core.layer_state) | verified |
| Cells layers with level of detail | yes | yes (same core) | not yet exercised with real segmentation |
| Spots layer | yes | yes (lazy, same core) | API path only |
| Colour by gene / metadata, categories | yes | yes | verified (gene); categories API-tested |
| Expression values (counts / normalized / log-normalized) | yes | yes | verified |
| Colour, opacity, show not sequenced, fixed bin resolution | yes | yes | UI wired; not all automated |
| Hover barcode (fetched on demand, cached) | yes | yes | verified (Playwright) |
| Map selection (box / lasso) to native barcodes | yes | yes | API-tested; drives later analyses |
| Thumbnail with viewport box, copy viewport coordinates | yes | yes | rendered (Playwright) |
| Ruler | yes | yes (same Bokeh tool) | in the figure toolbar |
| Image opacity slider | yes | yes (browser-side) | verified (Playwright) |
| Clustering: Leiden, k-means, DBSCAN, HDBSCAN, Gaussian mixture, agglomerative, with their options | yes | yes (background job, live progress) | verified (Leiden in Playwright on 16 µm mouse brain; API tests) |
| Clustering on the displayed bin aggregation level | yes | yes (shared core.clustering) | verified |
| Cluster list in assigned colours, sizes, marker genes | yes | yes | verified |
| Show / hide clusters, select all / none | yes | yes (map filter) | verified |
| Recolour a cluster | via DE/legends | yes (API) | API-tested; no picker in the list yet |
| Groups from gene/metadata conditions (AND/OR, validation hints) | yes | yes | verified (Playwright + API) |
| Show / hide, recolour, delete groups; colour the map by groups | yes | yes | API-tested; create verified in browser |
| ROI by bounding box or lasso drawn on the map (own tools: drawing selects nothing) | yes | yes | verified (Playwright: box) + API (lasso) |
| ROI by seed + radius: click on map, barcode / bin index (suggestions), all bins of a group (full or minus the group) | yes | yes | verified (Playwright: map click) + API (barcode, group) |
| Radius in µm when the image is calibrated, else full-resolution px | yes | yes | API-tested |
| Draft preview, save / cancel; one drawing at a time across layers | yes | yes | verified (Playwright) |
| Group-seed previews debounced, superseded previews cancelled, save reuses the preview | yes | yes | API-tested |
| Exact membership on native geometries; large outlines simplified for display | yes (raster overview) | yes (simplified vector) | API-tested |
| ROI list: count, seeds, focus, select, unselect, create group, delete (with its group), outline toggle | yes | yes | verified (Playwright) + API |
| DE: groups and/or clusters, vs rest or between selected, resolution (displayed / native / aggregated) | yes | yes (background job) | verified (Playwright: group vs rest) + API (between groups) |
| DE table below the map: set colours, L2FC / p-adj colouring, search and sort over every gene | yes (200 rows client-side) | yes (server pages of 200, "Show more") | verified (Playwright) + API |
| DE plots (violins for bins / cells, dots for spots) | yes | yes (same Bokeh figures) | verified (Playwright) |
| DE results cleared when a compared group is deleted or the layer is re-clustered | yes | yes (a running DE is discarded too) | verified (Playwright) + API |
| Deconvolution methods from plugin manifests (NNLS, RCTD, MuSiC, cell2location, CIBERSORTx) with their parameter forms, incl. reference-column parameters | yes | yes | NNLS verified (API, real run); others need their containers |
| Reference cell-type column; reference summary and runner shown | yes | yes | API-tested |
| Live availability hint (no reference, runner missing, invalid parameters) | yes | yes (debounced check) | verified (Playwright: no reference) + API |
| Run as a background job with progress; separate from clustering/DE; cancelled when the session ends | yes | yes | API-tested (NNLS, local runner) |
| Results cached under the dataset per input resolution and layer (same folders as the Bokeh UI) | yes | yes (shared core.deconvolution) | API-tested |
| Saved runs, readout (fractions / abundance), display a cell type (first one shown after a run) | yes | yes | API-tested; browser: see e2e/deconvolution.spec.ts (Human CRC) |
| Metadata comparison: groups and/or clusters, or each against The rest; numeric metadata or a deconvolution readout; mean / median | yes | yes | verified (Playwright: numeric) + API (fractions, The rest) |
| Latest deconvolution run preselected as the metadata set | yes | yes | API-tested |
| Comparison cleared when a compared set is deleted, the layer is re-clustered, or (with The rest) any group changes | yes | yes | verified (Playwright) + API |
| Summary: dataset table, QC figures, marker genes after clustering | yes | yes | verified (API + Playwright, 16 µm mouse brain) |
| Embeddings (PCA, UMAP, t-SNE) coloured by a column or the clusters' colours, linked selection with the map and ROIs | yes | yes | verified (API + Playwright, 16 µm mouse brain) |
| Sample comparison (spots vs bins): metrics table, local view, show on map, linked selection | yes | yes | verified (API + Playwright, 16 µm mouse brain) |
| Download a layer's data (.h5ad with clusters, groups, deconvolution results; private temporary file removed after sending) | yes (non-lazy layers only) | yes (every layer) | verified (API + Playwright) |
| Parameters tab: fill opacity, line width, line colour, line opacity, fixed or from data (categories for colours, genes or numeric metadata scaled to a range), optional hover column | yes | yes (map rebuilt at the current view) | verified (API: all attributes; Playwright: line width) |
