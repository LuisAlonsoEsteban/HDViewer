
from bokeh.plotting import curdoc
from bokeh.models import Div

from .panels.startup import create_startup


doc = curdoc()
create_startup(doc)

placeholders = {
    "metadata": Div(text=""),
    "spots_widgets": Div(text=""),
    "bins_widgets": Div(text=""),
    "cells_widgets": Div(text=""),
    "plot_column": Div(text=""),
    "thumbnail_plot": Div(text=""),
    "spots_params": Div(text=""),
    "bins_params": Div(text=""),
    "cells_params": Div(text=""),
}

# Add them to the document with names
for name, comp in placeholders.items():
    comp.name = name
    doc.add_root(comp)
