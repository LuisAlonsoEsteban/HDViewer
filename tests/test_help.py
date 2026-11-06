"""Help texts (GET /api/help) and the tutorial page they link to."""
import re
import unittest
from html.parser import HTMLParser
from pathlib import Path

from fastapi.testclient import TestClient

from hd_viewer.api.main import create_app
from hd_viewer.api.sessions import SessionRegistry
from hd_viewer.help import BIOLOGICAL_USE_CASE_URL, HELP, TUTORIAL_ANCHORS, TUTORIAL_URL, help_topics

ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = ROOT / "hd_viewer" / "static" / "tutorial.html"
FRONTEND = ROOT / "frontend" / "src"


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


def frontend_topics():
    """Help topics the React interface shows ("?" buttons): {file name: [topic]}."""
    used = {}
    for path in FRONTEND.rglob("*.tsx"):
        text = path.read_text()
        topics = re.findall(r'topic="([a-z_]+)"', text)
        for group in re.findall(r"topics=\{\[([^\]]*)\]\}", text):
            topics += re.findall(r'"([a-z_]+)"', group)
        if topics:
            used[path.name] = topics
    return used


class TutorialPageTests(unittest.TestCase):
    def setUp(self):
        self.page = Ids()
        self.page.feed(TUTORIAL.read_text())

    def test_served_next_to_the_app(self):
        client = TestClient(create_app(SessionRegistry(ttl=3600), sweep_seconds=3600))
        for url, name in ((TUTORIAL_URL, TUTORIAL.name), (BIOLOGICAL_USE_CASE_URL, "biological-use-case.html")):
            self.assertFalse(url.startswith("/"))                # relative: works under any URL prefix
            self.assertEqual(url.rsplit("/", 1)[1], name)
            response = client.get("/" + url)
            self.assertEqual(response.status_code, 200, url)
            self.assertEqual(response.content, (ROOT / "hd_viewer" / "static" / name).read_bytes())

    def test_two_tabs_and_no_broken_internal_links(self):
        self.assertTrue({"load", "app", "tab-load", "tab-app"} <= self.page.ids)
        self.assertEqual(self.page.links - self.page.ids, set())

    def test_every_help_topic_links_to_an_existing_section(self):
        self.assertEqual(set(TUTORIAL_ANCHORS), set(HELP))
        topics = help_topics()
        for topic, anchor in TUTORIAL_ANCHORS.items():
            self.assertIn(anchor, self.page.ids, topic)
            self.assertEqual(topics[topic]["tutorial"], f"{TUTORIAL_URL}#{anchor}")


class HelpTests(unittest.TestCase):
    def test_help_endpoint_serves_every_topic(self):
        client = TestClient(create_app(SessionRegistry(ttl=3600), sweep_seconds=3600))
        served = client.get("/api/help").json()
        self.assertEqual(served, help_topics())
        self.assertEqual(set(served), set(HELP))
        for topic, entry in served.items():
            self.assertEqual(entry["title"], HELP[topic][0])
            self.assertTrue(entry["html"].strip(), topic)
            self.assertNotIn("\n", entry["html"])

    def test_loading_screen_explains_its_sections(self):
        self.assertEqual(sorted(HELP[t][0] for t in frontend_topics()["Loading.tsx"]),
                         sorted(["Your data", "What to compute", "Loading progress"]))

    def test_every_layer_tab_explains_itself_with_existing_topics(self):
        used = frontend_topics()
        every = [topic for topics in used.values() for topic in topics]
        self.assertEqual(set(every) - set(HELP), set())          # no button without a text
        self.assertEqual(set(HELP) - set(every), set())          # no text without a button
        self.assertIn("bin_resolution", used["BinInfoOverlay.tsx"])
        workspace = (FRONTEND / "screens" / "Workspace.tsx").read_text()
        panels = [panel for panel in re.findall(r"<Tabs\.Panel\b[^>]*>(.*)", workspace) if "layer={layer}" in panel]
        self.assertGreaterEqual(len(panels), 8)                  # the per-layer tabs
        for panel in panels:
            self.assertIn("<PanelHelp topics=", panel)


if __name__ == "__main__":
    unittest.main()
