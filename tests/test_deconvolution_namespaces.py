"""Human-readable cache namespaces follow native inputs, not display zoom."""
import unittest
from pathlib import Path
from types import SimpleNamespace

from hd_viewer.core.deconvolution import cache_namespace
from hd_viewer.startup.sample import Cells, Spots
from tests.test_deconvolution import DeconvolutionApiTestCase, reference
from tests.test_param_form import make_layer


class DeconvolutionNamespaceTests(unittest.TestCase):
    def test_native_bin_size_separates_cache_and_display_zoom_does_not(self):
        view = make_layer()
        for size, expected in ((2, 'square_002um/bins'), (8, 'square_008um/bins'), (16, 'square_016um/bins')):
            view.sample.bin_size = size
            self.assertEqual(cache_namespace(view.sample, view.kind, view.id), expected)
            view.set_aggregation(4)
            self.assertEqual(cache_namespace(view.sample, view.kind, view.id), expected)
            view.set_aggregation(1)

    def test_namespace_encodes_layer_identifiers_as_single_directories(self):
        sample = make_layer().sample
        namespace = Path(cache_namespace(sample, 'cells', '../../other/layer', 2))
        self.assertEqual(namespace.parts[:2], ('square_002um', 'cells'))
        self.assertEqual(len(namespace.parts), 3)
        self.assertNotIn('..', namespace.parts)
        self.assertEqual(cache_namespace(sample, 'cells', 'x', 0.5), 'square_0p5um/cells/x')
        self.assertEqual(cache_namespace(SimpleNamespace(), 'spots', 'spots'), 'native/spots')


class LayerNamespaceTests(DeconvolutionApiTestCase):
    def test_cells_and_spots_keep_input_resolution_and_distinct_layer_names(self):
        bins = make_layer().sample
        samples = {'cells_method_a': ('cells', Cells(bins.gdf.copy(), bins.adata.copy()), 'Cells · a'),
                   'cells_method_b': ('cells', Cells(bins.gdf.copy(), bins.adata.copy()), 'Cells · b'),
                   'spots': ('spots', Spots(bins.gdf.copy(), bins.adata.copy()), 'Spots')}
        base = Path(self.folder.name) / 'dataset'
        session = self.install(samples, reference(), base, dataset=SimpleNamespace(bins_size=8))
        for layer, suffix in (('cells_method_a', 'cells/cells_method_a'), ('cells_method_b', 'cells/cells_method_b'),
                              ('spots', 'spots')):
            with self.subTest(layer=layer):
                self.run_method(layer)
                uns = session.views[layer].sample.adata.uns['deconvolution']
                self.assertEqual(uns['cache_namespace'], f'square_008um/{suffix}')
                self.assertTrue(Path(uns['log_path']).is_relative_to(base / 'deconvolution' / 'square_008um' / suffix))


if __name__ == '__main__':
    unittest.main()
