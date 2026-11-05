"""Shared project identity and styles for Bokeh's individual shadow roots."""
from pathlib import Path

from bokeh.model import Model
from bokeh.document.events import ModelChangedEvent, RootAddedEvent
from bokeh.models import ImportedStyleSheet, UIElement
from bokeh.themes import Theme


BRANDING = {
    'name': 'VisiumVisualizer',
    'repository': 'https://github.com/LuisAlonsoEsteban/VisiumVisualizer',
    'issues': 'https://github.com/LuisAlonsoEsteban/VisiumVisualizer/issues',
    'author': 'Luis Alonso',
    'email': 'laesteban@unav.es',
    'affiliation': 'Computational Biology · Tecnun · University of Navarra',
    'affiliation_url': 'https://en.unav.edu/web/biomedical-engineering-centre/areas-of-research/biologia-computacional-e-ia',
    'citation': 'Alonso, L., et al. VisiumVisualizer [Computer software]. https://github.com/LuisAlonsoEsteban/VisiumVisualizer',
    'data_paper': 'https://doi.org/10.1038/s41588-025-02193-3',
}


def configure_document(doc):
    doc.title = BRANDING['name']
    doc.template_variables['branding'] = dict(BRANDING)
    if getattr(doc, '_visium_stylesheet', None) is not None:
        return
    # A fresh model per document: Bokeh models cannot belong to two sessions.
    sheet = doc._visium_stylesheet = ImportedStyleSheet(url='/hd_viewer/static/css/widgets.css')
    doc.theme = Theme(filename=str(Path(__file__).with_name('theme.yaml')))

    def style_models(value):
        models = set()
        values = value if isinstance(value, (list, tuple)) else [value]
        for model in values:
            if isinstance(model, Model):
                models.update(model.references())
        for model in models:
            if isinstance(model, UIElement) and not any(
                    isinstance(existing, ImportedStyleSheet) and existing.url == sheet.url
                    for existing in model.stylesheets):
                model.stylesheets = [*model.stylesheets, sheet]

    def changed(event):
        if isinstance(event, RootAddedEvent):
            style_models(event.model)
        elif isinstance(event, ModelChangedEvent) and event.attr in ('children', 'tabs'):
            # Only structural UI changes. Streaming data never walks the UI tree.
            # Styled on the next tick: this callback runs before the server session
            # has sent `event`, so changing the new models now would send updates
            # that reference models the browser does not know yet ("reference ...
            # isn't known"), and the browser would drop them.
            new = event.new
            doc.add_next_tick_callback(lambda: style_models(new))

    doc.on_change(changed)
    for root in doc.roots:
        style_models(root)
