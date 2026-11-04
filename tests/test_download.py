"""Exports are written on demand and contain the latest layer data."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
from bokeh.document import Document
from bokeh.events import ButtonClick

from tests.test_roi_manager import ImmediateDoc, make_layer


class DownloadTests(unittest.TestCase):
    def test_widget_builds_do_not_export_and_clicks_capture_current_data(self):
        with tempfile.TemporaryDirectory() as root:
            folder = Path(root)/'exports'
            with patch('hd_viewer.board.components.layer._layer.DOWNLOADABLE_DATA_PATH', folder):
                component = make_layer('spots')
                self.assertFalse(folder.exists())
                component.get_widgets(ImmediateDoc())
                self.assertFalse(folder.exists())
                button = component.download_button
                self.assertEqual(button.tags, [])
                component.sample.adata.obs['latest_result'] = 0.25
                button._trigger_event(ButtonClick(button))
                first = component.adata_path
                self.assertEqual(len(list(folder.glob('*.h5ad'))), 1)
                restored = ad.read_h5ad(first)
                np.testing.assert_allclose(restored.obs.latest_result, .25)
                self.assertEqual(button.tags[0]['url'], component.static_adata_url)
                # Exports live in a public folder: names must be unique and unguessable.
                self.assertRegex(first.name, r'^spots-[0-9a-f]{32}\.h5ad$')
                self.assertEqual(button.tags[0]['filename'], 'spots.h5ad')
                component.get_widgets(ImmediateDoc())
                self.assertEqual(len(list(folder.glob('*.h5ad'))), 1)
                component.sample.adata.obs['latest_result'] = 0.75
                component.download_button._trigger_event(ButtonClick(component.download_button))
                self.assertEqual(len(list(folder.glob('*.h5ad'))), 2)
                self.assertNotEqual(first, component.adata_path)
                np.testing.assert_allclose(ad.read_h5ad(component.adata_path).obs.latest_result, .75)
                doc = Document()
                doc.add_root(component.download_button)
                doc.to_json()

    def test_failed_export_is_removed_and_can_be_retried(self):
        with tempfile.TemporaryDirectory() as root, \
             patch('hd_viewer.board.components.layer._layer.DOWNLOADABLE_DATA_PATH', Path(root)):
            component = make_layer('spots')
            button = component.download_button
            with patch.object(component.sample.adata, 'write_h5ad', side_effect=OSError('test failure')):
                with self.assertRaises(OSError):
                    button._trigger_event(ButtonClick(button))
            self.assertFalse(list(Path(root).glob('*.h5ad')))
            self.assertEqual(button.tags, [])
            self.assertFalse(button.disabled)
            button._trigger_event(ButtonClick(button))
            self.assertEqual(len(list(Path(root).glob('*.h5ad'))), 1)


if __name__ == '__main__':
    unittest.main()
