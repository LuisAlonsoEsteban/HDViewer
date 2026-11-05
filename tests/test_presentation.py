"""Presentation keeps dataset controls reachable and Bokeh sessions isolated."""
from html.parser import HTMLParser
from pathlib import Path
from types import SimpleNamespace
import unittest
from urllib.parse import urlsplit

from bokeh.document import Document
from bokeh.layouts import column
from bokeh.models import Button, ColumnDataSource, ImportedStyleSheet, InlineStyleSheet
from jinja2 import Environment, FileSystemLoader, Template
import yaml

from hd_viewer.presentation import BRANDING, configure_document


ROOT = Path(__file__).resolve().parents[1]


class PageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.elements = []

    def handle_starttag(self, tag, attributes):
        self.elements.append((tag, dict(attributes)))


class PresentationTests(unittest.TestCase):
    def render(self):
        env = Environment(loader=FileSystemLoader(ROOT/'hd_viewer'/'templates'))
        base = Template('{% block preamble %}{% endblock %}{% block contents %}{% endblock %}')
        names = ('startup', 'metadata', 'embedding_panel', 'spots_widgets', 'bins_widgets',
                 'cells_widgets', 'comparison_panel', 'plot_column', 'de_results', 'thumbnail_plot')
        html = env.get_template('index.html').render(base=base, branding=BRANDING,
                    layout_settings={}, roots=SimpleNamespace(**{name:name for name in names}),
                    embed=lambda name:f'<div data-bokeh-root="{name}"></div>')
        parser = PageParser()
        parser.feed(html)
        return html, parser.elements

    def test_landing_keeps_all_analysis_roots_and_setup_in_one_session(self):
        html, elements = self.render()
        ids = {attrs['id']:attrs for _, attrs in elements if 'id' in attrs}
        self.assertNotIn('hidden', ids['landing-container'])
        self.assertIn('hidden', ids['startup-container'])
        self.assertIn('hidden', ids['main-container'])
        self.assertEqual(len([attrs for _, attrs in elements if 'data-bokeh-root' in attrs]), 10)
        self.assertIn('software-citation', ids)
        self.assertIn('Cima Universidad de Navarra', html)
        self.assertIn('Tecnun Universidad de Navarra', html)
        self.assertIn(BRANDING['repository'], html)
        self.assertIn('mailto:'+BRANDING['email'], html)

    def test_local_asset_links_exist_and_external_links_are_safe(self):
        _, elements = self.render()
        for tag, attrs in elements:
            target = attrs.get('src', attrs.get('href', ''))
            if target.startswith('/hd_viewer/static/'):
                self.assertTrue((ROOT/'hd_viewer'/'static'/urlsplit(target).path.split('/static/',1)[1]).is_file(), target)
            if attrs.get('target') == '_blank':
                self.assertIn('noopener', attrs.get('rel', ''), target)
            if tag == 'img':
                self.assertTrue(attrs.get('alt'))

    def test_shared_styles_are_owned_by_each_document(self):
        first, second = Document(), Document()
        for doc in (first, second):
            configure_document(doc)
            doc.add_root(column(Button(label='Test')))
            self.assertEqual(doc.title, BRANDING['name'])
            doc.to_json()
        self.assertIsNot(first._visium_stylesheet, second._visium_stylesheet)
        self.assertEqual(first.roots[0].children[0].stylesheets, [first._visium_stylesheet])

    def test_dynamic_widgets_keep_custom_styles_and_receive_shared_styles_once(self):
        from tests.test_patch_ordering import ServerLikeDoc
        doc = ServerLikeDoc()
        configure_document(doc)
        layout = column()
        doc.add_root(layout)
        original = InlineStyleSheet(css=':host { margin: 1px; }')
        button = Button(stylesheets=[original])
        layout.children = [button]
        # Styled on the next tick, once the new widgets have reached the browser
        self.assertEqual(button.stylesheets, [original])
        doc.pump(.1)
        configure_document(doc)
        self.assertIs(button.stylesheets[0], original)
        self.assertEqual(sum(isinstance(sheet, ImportedStyleSheet) for sheet in button.stylesheets), 1)
        layout.children = [button, Button()]
        doc.pump(.1)
        self.assertEqual(len(button.stylesheets), 2)
        self.assertEqual(layout.children[1].stylesheets, [doc._visium_stylesheet])

    def test_downloadable_citation_matches_repository_metadata(self):
        original = (ROOT/'CITATION.cff').read_text()
        self.assertEqual(original, (ROOT/'hd_viewer'/'static'/'citation.cff').read_text())
        citation = yaml.safe_load(original)
        self.assertEqual(citation['repository-code'], BRANDING['repository'])
        self.assertEqual(citation['title'], BRANDING['name'])
        self.assertNotIn('doi', citation)

    def test_session_callbacks_and_streamed_data_do_not_restyle_the_tree(self):
        doc = Document()
        configure_document(doc)
        root = column(Button())
        doc.add_root(root)
        callback = doc.add_periodic_callback(lambda: None, 100)
        doc.remove_periodic_callback(callback)
        source = ColumnDataSource({'x': [1.]})
        doc.add_root(source)
        from unittest.mock import patch
        with patch.object(type(root), 'references', side_effect=AssertionError('Data update walked the UI tree')):
            source.stream({'x': [2.]})


if __name__ == '__main__':
    unittest.main()
