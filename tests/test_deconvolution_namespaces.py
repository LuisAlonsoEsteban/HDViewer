"""Human-readable cache namespaces follow native inputs, not display zoom."""
import unittest
from pathlib import Path
from bokeh.document import Document

import hd_viewer.startup
from hd_viewer.board.app import App
from tests.test_roi_manager import make_layer


class DeconvolutionNamespaceTests(unittest.TestCase):
    def test_native_bin_size_separates_cache_and_display_zoom_does_not(self):
        layer = make_layer()
        manager = layer.deconvolution_manager
        for size, expected in ((2, 'square_002um/bins'), (8, 'square_008um/bins'), (16, 'square_016um/bins')):
            layer.sample.bin_size = size
            self.assertEqual(manager.cache_namespace, expected)
            layer.manager.set_aggregation_factor(4)
            self.assertEqual(manager.cache_namespace, expected)

    def test_cells_and_spots_keep_input_resolution_and_distinct_layer_names(self):
        app = App(Document())
        app.deconvolution_bin_size_um = 8
        app.deconvolution_work_dir = Path('/tmp/example-deconvolution')
        for kind, suffix in (('method_a', 'cells/cells_method_a'), ('method_b', 'cells/cells_method_b'),
                             ('spots', 'spots')):
            layer = make_layer(kind)
            app.add_geometry(layer)
            self.assertEqual(layer.deconvolution_manager.cache_namespace, f'square_008um/{suffix}')
            self.assertEqual(layer.deconvolution_manager.work_dir, app.deconvolution_work_dir)

    def test_namespace_encodes_layer_identifiers_as_single_directories(self):
        layer = make_layer('method_a')
        layer.deconvolution_manager.source_bin_size_um = 2
        layer.key = '../../other/layer'
        namespace = Path(layer.deconvolution_manager.cache_namespace)
        self.assertEqual(namespace.parts[:2], ('square_002um', 'cells'))
        self.assertEqual(len(namespace.parts), 3)
        self.assertNotIn('..', namespace.parts)


if __name__ == '__main__':
    unittest.main()
