"""Check shared configuration boundaries and non-default settings."""
import inspect
import json
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from jinja2 import Environment

from hd_viewer import settings
from hd_viewer.settings.config import sample_kwargs, image_kwargs, deconvolution_kwargs
from hd_viewer.board.components.image import ImageComponent
from hd_viewer.board.layout.layer.parameter import attr_spec
from hd_viewer.board.source.lazyloader import LazyLoader
from hd_viewer.startup import StartupController
from hd_viewer.startup.data.image.image import ImageProcessor
from hd_viewer.startup.sample import Bins, Cells
from hd_viewer.startup.sample._sample import Sample
from hd_viewer.startup.data.image.segmentation import plugins
from hd_viewer.board.operations.layer.deconvolution import deconvolve
from tests.test_roi_manager import ImmediateDoc, make_layer


class SettingsTests(unittest.TestCase):
    def test_shared_settings_and_compatibility_aliases(self):
        self.assertEqual(ImageComponent.TILE_MEM_THRESHOLD, settings.TILE_MEM_THRESHOLD)
        self.assertEqual(LazyLoader.OFFSET, settings.OFFSET)
        self.assertEqual(Bins.AGGREGATION_FACTORS, settings.AGGREGATION_FACTORS)
        self.assertEqual(Cells.LOD_MAX_LEVEL, settings.LOD_MAX_LEVEL)
        self.assertEqual(plugins.SDK_DIR, settings.SDK_DIR)
        self.assertEqual(image_kwargs()['tile_size'], settings.TILE_SIZE)
        self.assertEqual(deconvolution_kwargs()['batch_size'], settings.DECONV_BATCH_SIZE)
        self.assertEqual(settings.DOWNLOADABLE_DATA_PATH, Path(settings.__file__).parents[1]/'static'/'data')

    def test_custom_roi_radius_and_parameter_specs(self):
        with patch.object(settings, 'ROI_DEFAULT_RADIUS', 73.), \
             patch.object(settings, 'PARAMETER_WIDTH_SPEC', dict(low=0,high=12,value=2.5,not_sequenced=.2)):
            component = make_layer()
            self.assertEqual(component.roi_manager.radius.value, 73.)
            spec = attr_spec('line_width')
            self.assertEqual(spec['value'], 2.5)
            spec['value'] = 8
            self.assertEqual(settings.PARAMETER_WIDTH_SPEC['value'], 2.5)

    def test_configured_cluster_column_is_used_by_scanpy_and_markers(self):
        data = ad.AnnData(np.ones((4,2)),obs=pd.DataFrame(index=['a','b','c','d']),var=pd.DataFrame(index=['g','h']))
        def leiden(target, *, resolution, key_added, random_state):
            target.obs[key_added] = ['0','0','1','1']
        with patch.object(settings, 'CLUSTER_COL', 'custom_leiden'):
            sample = Sample(pd.DataFrame({'index':list(data.obs_names)}), data, **sample_kwargs())
        sample.adata_processed = data.copy()
        with patch.object(sample, '_ensure_neighbors', return_value=sample.adata_processed), \
             patch('scanpy.tl.leiden', side_effect=leiden), \
             patch('scanpy.tl.rank_genes_groups') as markers:
            sample.run_clustering(key_added='custom_leiden')
            self.assertEqual(sample.adata.obs.custom_leiden.tolist(), ['0','0','1','1'])
            sample.compute_marker_genes()
            self.assertEqual(markers.call_args.kwargs['groupby'], 'custom_leiden')

    def test_layout_settings_reach_template(self):
        layout = {**settings.PANEL_SPLITS,'map':dict(sizes=[60,40],minSize=[150,70],gutterSize=9)}
        with patch.object(settings,'PANEL_SPLITS',layout):
            controller = StartupController(ImmediateDoc())
            rendered = Environment().from_string('{{ layout_settings|tojson }}').render(**controller.doc.template_variables)
            self.assertEqual(json.loads(rendered)['map']['sizes'], [60,40])


if __name__ == '__main__':
    unittest.main()
