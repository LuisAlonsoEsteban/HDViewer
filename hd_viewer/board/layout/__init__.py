
from bokeh.layouts import row, column
from bokeh.models import Div

class Panel:

    def __init__(
            self, 
            orient='col', 
            name='', 
            sizing_mode='stretch_both'):
        
        if orient == 'row':
            self.layout = row(name=name)
        elif orient == 'col':
            self.layout = column(name=name)

        self.set_placeholder()

    def add_children(self, children):
        if isinstance(children, list):
            self.layout.children.extend(children)
        else:
            self.layout.children.append(children)
        
    def set_children(self, children):
        self.layout.children = children

    def get_placeholder(self):
        return [Div(
            text="Empty panel"
        )]

    def set_placeholder(self):
        self.set_children(self.get_placeholder())

    def get_layout(self):
        return self.layout


from .thumbnail import ThumbnailPanel
from .figure import FigurePanel
from .metadata import MetadataPanel
from .embedding import EmbeddingPanel
from .comparison import ComparisonPanel
