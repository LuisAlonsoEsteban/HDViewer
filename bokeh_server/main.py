
from bokeh.plotting import curdoc
from bokeh.layouts import column

from .panels.startup import create_startup


doc = curdoc()
create_startup(doc)

placeholders = {
    "metadata": column(sizing_mode='stretch_both'),
    "spots_widgets": column(sizing_mode='stretch_both'),
    "bins_widgets": column(sizing_mode='stretch_both'),
    "cells_widgets": column(sizing_mode='stretch_both'),
    "plot_column": column(sizing_mode='stretch_both'),
    "thumbnail_plot": column(sizing_mode='stretch_both'),
    "spots_params": column(sizing_mode='stretch_both'),
    "bins_params": column(sizing_mode='stretch_both'),
    "cells_params": column(sizing_mode='stretch_both'),
}

# Add them to the document with names
for name, comp in placeholders.items():
    comp.name = name
    doc.add_root(comp)
