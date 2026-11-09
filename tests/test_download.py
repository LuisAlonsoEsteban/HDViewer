"""Layer exports (core.export, behind GET /api/layers/{id}/export) are written on
demand, hold the latest layer data and never leave partial files behind."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np

from hd_viewer.core.export import export_layer
from tests.test_param_form import make_layer


class ExportTests(unittest.TestCase):
    def test_each_export_captures_current_data_in_a_new_unguessable_file(self):
        with tempfile.TemporaryDirectory() as root:
            folder = Path(root) / 'exports'
            view = make_layer()
            self.assertFalse(folder.exists())               # nothing is written until asked
            view.sample.adata.obs['latest_result'] = 0.25
            first = export_layer(view.sample, folder)
            self.assertEqual(first.parent, folder)
            self.assertRegex(first.name, r'^[0-9a-f]{32}\.h5ad$')
            np.testing.assert_allclose(ad.read_h5ad(first).obs.latest_result, .25)
            view.sample.adata.obs['latest_result'] = 0.75
            second = export_layer(view.sample, folder)
            self.assertNotEqual(first, second)
            self.assertEqual(len(list(folder.glob('*.h5ad'))), 2)
            np.testing.assert_allclose(ad.read_h5ad(second).obs.latest_result, .75)
            np.testing.assert_allclose(ad.read_h5ad(first).obs.latest_result, .25)

    def test_failed_export_is_removed_and_can_be_retried(self):
        with tempfile.TemporaryDirectory() as root:
            view = make_layer()
            with patch.object(ad.AnnData, 'write_h5ad', side_effect=OSError('test failure')):
                with self.assertRaises(OSError):
                    export_layer(view.sample, root)
            self.assertFalse(list(Path(root).glob('*.h5ad')))
            path = export_layer(view.sample, root)
            self.assertEqual(list(Path(root).glob('*.h5ad')), [path])
            self.assertEqual(ad.read_h5ad(path).shape, view.sample.adata.shape)

    def test_unstorable_metadata_falls_back_without_losing_the_data(self):
        with tempfile.TemporaryDirectory() as root:
            view = make_layer()
            view.sample.adata.uns['unstorable'] = object()
            view.sample.adata.obs['mixed'] = ['a', 1, None, 2.5, 'b', 'c']
            restored = ad.read_h5ad(export_layer(view.sample, root))
            self.assertNotIn('unstorable', restored.uns)
            self.assertEqual(list(restored.obs['mixed'].astype(object).where(restored.obs['mixed'].notna(), None)),
                             ['a', '1', None, '2.5', 'b', 'c'])
            self.assertIn('unstorable', view.sample.adata.uns)        # the layer itself is unchanged
            self.assertEqual(len(list(Path(root).glob('*.h5ad'))), 1)


if __name__ == '__main__':
    unittest.main()
