"""Check shared configuration boundaries and non-default settings."""
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

from hd_viewer import settings
from hd_viewer.core.clustering import ClusteringState
from hd_viewer.core.parameters import attr_spec
from hd_viewer.settings.config import deconvolution_kwargs, image_kwargs, sample_kwargs
from hd_viewer.startup.data.image.segmentation import plugins
from hd_viewer.startup.sample import Bins, Cells
from hd_viewer.startup.sample._sample import Sample
from tests.fixtures import make_view, show


class SettingsTests(unittest.TestCase):
    def test_shared_settings_and_compatibility_aliases(self):
        self.assertEqual(Bins.AGGREGATION_FACTORS, settings.AGGREGATION_FACTORS)
        self.assertEqual(Cells.LOD_MAX_LEVEL, settings.LOD_MAX_LEVEL)
        self.assertEqual(plugins.SDK_DIR, settings.SDK_DIR)
        self.assertEqual(image_kwargs()['tile_size'], settings.TILE_SIZE)
        self.assertEqual(deconvolution_kwargs()['batch_size'], settings.DECONV_BATCH_SIZE)

    def test_viewport_settings_reach_the_layer_view(self):
        for offset in (settings.OFFSET, 1.5):
            with self.subTest(offset=offset), patch.object(settings, 'OFFSET', offset):
                view = make_view()
                show(view, 0, 2, 0, 2)
                margin = 2 * (1 + offset) / 2
                self.assertEqual(view.rows.bbox, [1 - margin, 1 + margin, 1 - margin, 1 + margin])

    def test_tile_memory_threshold_chooses_single_image(self):
        from hd_viewer.core.pipeline import make_image_processor
        from tests.test_spatial_alignment import fixture
        _, transform, _, _ = fixture()
        dataset = type('Dataset', (), dict(img=np.zeros((140, 180, 3), np.uint8), img_path=__file__))()
        with tempfile.TemporaryDirectory() as root:
            for threshold, single in ((0, False), (10 ** 12, True)):
                with patch.object(settings, 'TILE_MEM_THRESHOLD', threshold):
                    self.assertEqual(make_image_processor(dataset, transform, root)[1], single)

    def test_custom_roi_radius_and_parameter_specs(self):
        from fastapi import FastAPI
        from fastapi.testclient import TestClient
        from hd_viewer.api.routers import config
        app = FastAPI()
        app.include_router(config.router, prefix='/api')
        with patch.object(settings, 'ROI_DEFAULT_RADIUS', 73.), \
             patch.object(settings, 'PARAMETER_WIDTH_SPEC', dict(low=0,high=12,value=2.5,not_sequenced=.2)):
            self.assertEqual(TestClient(app).get('/api/config').json()['roi']['default_radius'], 73.)
            spec = attr_spec('line_width')
            self.assertEqual(spec['value'], 2.5)
            spec['value'] = 8
            self.assertEqual(settings.PARAMETER_WIDTH_SPEC['value'], 2.5)

    def test_custom_cluster_column_is_used_by_scanpy_and_markers(self):
        data = ad.AnnData(np.ones((4,2)),obs=pd.DataFrame(index=['a','b','c','d']),var=pd.DataFrame(index=['g','h']))
        def leiden(target, *, resolution, key_added, random_state):
            target.obs[key_added] = ['0','0','1','1']
        sample = Sample(pd.DataFrame({'index':list(data.obs_names)}), data, **sample_kwargs())
        clustering = ClusteringState(sample)
        clustering.adata_processed = data.copy()
        with patch.object(clustering, '_ensure_neighbors', return_value=clustering.adata_processed), \
             patch('scanpy.tl.leiden', side_effect=leiden), \
             patch('scanpy.tl.rank_genes_groups') as markers:
            clustering.run_clustering(key_added='custom_leiden')
            self.assertEqual(sample.adata.obs.custom_leiden.tolist(), ['0','0','1','1'])
            clustering.compute_markers()
            self.assertEqual(markers.call_args.kwargs['groupby'], 'custom_leiden')


if __name__ == '__main__':
    unittest.main()
