"""The load pipeline without any UI: progress per task and the samples it builds."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from PIL import Image

from hd_viewer.core.pipeline import LoadConfig, LoadPipeline
from hd_viewer.startup.sample import Bins
from tests.test_deconvolution import tenx_reference_bytes
from tests.test_spatial_alignment import fixture


def write_dataset(root):
    positions, _, _, _ = fixture()
    positions = positions.rename(index={'b0': 'a', 'b1': 'b', 'b6': 'c'})
    image = root / 'image.png'
    Image.fromarray(np.full((140, 180, 3), 180, dtype=np.uint8)).save(image)
    counts = root / 'counts.h5'
    counts.write_bytes(tenx_reference_bytes())
    positions_path = root / 'positions.parquet'
    positions.drop(columns='index').rename_axis('barcode').reset_index().to_parquet(positions_path)
    (root / 'scalefactors_json.json').write_text(json.dumps({'spot_diameter_fullres': 4}))
    return dict(image_path=str(image), expression_path=str(counts), positions_path=str(positions_path), bins_size=2)


class LoadPipelineTests(unittest.TestCase):
    def test_example_load_reports_every_task_and_builds_layers(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            config = LoadConfig.from_example(write_dataset(root), pseudospots=False, owner='session')
            # Three bins are too few for a 55 µm pseudospot: spots are covered elsewhere
            self.assertEqual(config.tasks, ['tiles', 'bins'])
            events = []
            loaded = LoadPipeline(config, progress=lambda task, pct, status: events.append((task, pct, status))).run()

            self.assertEqual({task for task, _, _ in events}, {'tiles', 'bins'})
            for task in config.tasks:
                self.assertIn((task, 100), [(t, pct) for t, pct, _ in events if pct == 100])
            self.assertEqual(set(loaded.samples), {'bins'})
            kind, bins, title = loaded.samples['bins']
            self.assertIsInstance(bins, Bins)
            self.assertIs(bins.spatial_transform, loaded.transform)
            self.assertEqual(bins.adata.n_obs, 3)
            self.assertTrue((loaded.image_processor.tiles_path / 'manifest.json').exists())

            # A second load reuses the cached tiles
            again = []
            LoadPipeline(config, progress=lambda task, pct, status: again.append((task, status))).run()
            self.assertIn(('tiles', 'Using cached tiles'), again)


if __name__ == '__main__':
    unittest.main()
