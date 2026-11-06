"""Project identity: the downloadable citation matches the repository metadata."""
from pathlib import Path
import unittest

import yaml

from hd_viewer.core.branding import BRANDING

ROOT = Path(__file__).resolve().parents[1]


class CitationTests(unittest.TestCase):
    def test_downloadable_citation_matches_repository_metadata(self):
        original = (ROOT/'CITATION.cff').read_text()
        self.assertEqual(original, (ROOT/'hd_viewer'/'static'/'citation.cff').read_text())
        citation = yaml.safe_load(original)
        self.assertEqual(citation['repository-code'], BRANDING['repository'])
        self.assertEqual(citation['title'], BRANDING['name'])
        self.assertNotIn('doi', citation)


if __name__ == '__main__':
    unittest.main()
