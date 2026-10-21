import re
import unittest
from html.parser import HTMLParser
from pathlib import Path

from bokeh.document import Document
from bokeh.models import HelpButton

from bokeh_server.help import HELP, TUTORIAL_ANCHORS, TUTORIAL_URL, help_button
from bokeh_server.panels.startup import StartupController
from tests.test_param_form import make_layer

TUTORIAL = Path(__file__).resolve().parents[1] / "bokeh_server" / "static" / "tutorial.html"


class Ids(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids, self.links = set(), set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("href", "").startswith("#"):
            self.links.add(attrs["href"][1:])


def topics_in(model):
    """Titles of the help buttons inside a Bokeh model."""
    return [re.search(r"<b style='font-size:13px'>(.*?)</b>", "".join(m.tooltip.content.html)).group(1)
            for m in model.references() if isinstance(m, HelpButton)]


class TutorialPageTests(unittest.TestCase):
    def setUp(self):
        self.page = Ids()
        self.page.feed(TUTORIAL.read_text())

    def test_served_next_to_the_app(self):
        self.assertTrue(TUTORIAL_URL.startswith("/bokeh_server/static/"))
        self.assertEqual(TUTORIAL_URL.rsplit("/", 1)[1], TUTORIAL.name)

    def test_two_tabs_and_no_broken_internal_links(self):
        self.assertTrue({"load", "app", "tab-load", "tab-app"} <= self.page.ids)
        self.assertEqual(self.page.links - self.page.ids, set())

    def test_every_help_button_links_to_an_existing_section(self):
        for topic, anchor in TUTORIAL_ANCHORS.items():
            self.assertIn(anchor, self.page.ids, topic)
            self.assertIn(f"{TUTORIAL_URL}#{anchor}", "".join(help_button(topic).tooltip.content.html))
        self.assertEqual(set(TUTORIAL_ANCHORS), set(HELP))


class HelpTests(unittest.TestCase):
    def test_startup_links_to_the_tutorial_and_explains_its_sections(self):
        startup = StartupController(Document())
        texts = [getattr(m, "text", "") for m in startup.layout.references()]
        self.assertTrue(any(TUTORIAL_URL in t and "Open the tutorial" in t for t in texts))
        self.assertFalse(any("<details>" in t for t in texts))          # no inline tutorial any more
        self.assertEqual(sorted(topics_in(startup.layout)),
                         sorted(["Your data", "What to compute", "Loading progress"]))

    def test_every_layer_tab_explains_itself(self):
        layer, _ = make_layer()
        tabs = layer.get_widgets(Document())
        for tab in tabs.tabs:
            self.assertEqual(len(topics_in(tab.child)), 1, tab.title)
        self.assertEqual(topics_in(layer.figure.elements[0]), ["Bin resolution"])


if __name__ == "__main__":
    unittest.main()
