
from bokeh.plotting import curdoc
from bokeh.layouts import column

from .panels.startup import create_startup


doc = curdoc()
create_startup(doc)

placeholders = {
    "metadata": column(),
    "embedding_panel": column(),
    "spots_widgets": column(),
    "bins_widgets": column(),
    "cells_widgets": column(),
    "plot_column": column(sizing_mode='stretch_width'),
    "thumbnail_plot": column(sizing_mode='stretch_both'),
}

# Add them to the document with names
for name, comp in placeholders.items():
    comp.name = name
    doc.add_root(comp)
