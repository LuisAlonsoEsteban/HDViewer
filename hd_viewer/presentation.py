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
    'author': 'Luis Alonso Esteban',
    'email': 'laesteban@unav.es',
    'affiliation': 'Computational Biology · Tecnun · University of Navarra',
    'affiliation_url': 'https://en.unav.edu/web/biomedical-engineering-centre/areas-of-research/biologia-computacional-e-ia',
    'citation': 'Alonso Esteban, L., et al. VisiumVisualizer [Computer software]. https://github.com/LuisAlonsoEsteban/VisiumVisualizer',
    'data_paper': None,
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
            style_models(event.new)

    doc.on_change(changed)
    for root in doc.roots:
        style_models(root)
