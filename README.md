# VisiumVisualizer

An interactive workspace for exploring Visium HD expression, cell populations,
and tissue neighbourhoods alongside the original tissue image.

Start the app from the repository root with `bokeh serve hd_viewer --show` (using
the environment where `requirements.txt` is installed). The landing page introduces
the project and opens the existing dataset-loading workflow. Returning to the
overview preserves the active workspace; its buttons then resume that session.

The landing page includes University of Navarra, Tecnun, and CIMA branding,
public project contact, source/issue links, and a repository-based software
citation. `CITATION.cff` and its downloadable copy at
`hd_viewer/static/citation.cff` should be kept in sync. The example's source-data
paper is cited separately from the software.

Presentation files:

- `hd_viewer/templates/index.html`: landing page and workspace shell.
- `hd_viewer/static/css/app.css`: responsive landing/setup/workspace layout.
- `hd_viewer/static/css/widgets.css`: shared styles inside Bokeh shadow roots.
- `hd_viewer/presentation.py`: project metadata and per-session style integration.
- `hd_viewer/static/assets/README.md`: image/logo sources and attribution.

## Deconvolution

The Deconvolution tab supports manifest-based NNLS, cell2location, MuSiC, and
CIBERSORTx plugins, with per-method controls and saved runs. See
[setup, method requirements, and plugin API](hd_viewer/board/operations/layer/deconvolution/api/README.md).

## Biological use-case walkthrough

The [`P1CRC walkthrough`](hd_viewer/static/biological-use-case.html) follows
Figure 4a–b of Oliveira et al. (2025), using `full_text.pdf` and the authors'
figure scripts. It covers 8 µm bins, Level2 RCTD labels, the 50 µm tumor periphery,
C1QC/COL1A1 maps, and composition comparisons. It explains the differences between
the viewer workflow and the published singlet-bin proportions.

Edit `hd_viewer/static/biological-use-case.html` to update the content. Give each
new section a unique `id` and add its link to the table of contents in that file.

For screenshots, create `hd_viewer/static/use-cases/` and reference images with
relative paths such as `use-cases/result.png`. Shared styling lives in
`hd_viewer/static/tutorial.css`.

The guide is linked from the loading screen, Summary panel, and tutorial. When the
app is running, open `/hd_viewer/static/biological-use-case.html` to preview it.
