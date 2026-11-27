"""Public-server behaviour: session creation and expiry, removal of orphaned
session folders, input validation and paths in messages; plus the faster
pooling and DE paths, checked against the straightforward versions."""
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from fastapi.testclient import TestClient

from hd_viewer import settings
from hd_viewer.analysis.de import DEAnalysis
from hd_viewer.api import sessions as sessions_module
from hd_viewer.api.main import create_app
from hd_viewer.api.routers import public_message
from hd_viewer.api.sessions import COOKIE, SessionRegistry
from hd_viewer.core.export import _without_server_paths
from hd_viewer.core.session import session_upload_dir
from hd_viewer.plugins.manifest import Param
from hd_viewer.startup.operations.pooler import Pooler


class ServerTestCase(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.uploads = Path(folder.name) / "uploads"
        self.uploads.mkdir()
        patcher = patch.object(settings, "UPLOAD_DIR", self.uploads)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.registry = SessionRegistry(ttl=3600, max_sessions=5)

    def client(self):
        client = TestClient(create_app(self.registry, sweep_seconds=3600, warm_up=False))
        client.__enter__()
        self.addCleanup(client.__exit__, None, None, None)
        return client


class SessionTests(ServerTestCase):
    def test_requests_without_a_session_do_not_create_one(self):
        client = self.client()
        response = client.post("/api/layers/bins/viewport", json=dict(seq=1, x0=0, x1=1, y0=0, y1=1))
        self.assertEqual(response.status_code, 409)
        self.assertEqual(client.get("/api/workspace").status_code, 409)
        self.assertEqual(len(self.registry), 0)
        self.assertNotIn(COOKIE, client.cookies)

    def test_new_session_cookie_is_set_on_error_responses(self):
        client = self.client()
        response = client.post("/api/uploads", json=dict(field="image", filename="x.exe", size=10))
        self.assertEqual(response.status_code, 422)
        self.assertIn(COOKIE, response.cookies)
        self.assertIsNotNone(self.registry.get(response.cookies[COOKIE]))
        # The browser keeps that session instead of creating one per request
        client.get("/api/load")
        self.assertEqual(len(self.registry), 1)

    def test_sessions_without_a_page_expire_early(self):
        client = self.client()
        client.get("/api/load")
        scripted = self.registry.get(client.cookies[COOKIE])
        page = self.registry.create()
        page.socket_seen = True
        working = self.registry.create()
        working.load = dict(state="running")
        for session in (scripted, page, working):
            session.last_seen = time.monotonic() - sessions_module.PRESOCKET_SECONDS - 1
        self.assertEqual(self.registry.sweep(), [scripted.sid])

    def test_closing_a_session_cancels_its_load(self):
        session = self.registry.create()
        session.load = dict(state="running", cancelled=False)
        self.registry.close(session.sid)
        self.assertTrue(session.load["cancelled"])


class OrphanFolderTests(ServerTestCase):
    def test_startup_and_sweeps_remove_folders_of_unknown_sessions(self):
        orphan = self.uploads / ("a" * 32)
        (orphan / "image").mkdir(parents=True)
        shared = self.uploads / "deconvolution"
        shared.mkdir()
        self.client()                                  # startup
        self.assertFalse(orphan.exists())
        self.assertTrue(shared.exists())

        live = self.registry.create()
        session_upload_dir(live.sid).mkdir()
        from hd_viewer.core.session import purge_orphan_folders
        orphan.mkdir()
        purge_orphan_folders(self.registry.active_folders())
        self.assertFalse(orphan.exists())
        self.assertTrue(session_upload_dir(live.sid).exists())


class ValidationTests(unittest.TestCase):
    def test_non_numeric_parameters_are_value_errors(self):
        param = Param(name="epochs", type="int", default=10, min=1, max=100)
        for value in ([1], {"a": 1}, True):
            with self.assertRaises(ValueError):
                param.coerce(value)
        with self.assertRaises(ValueError):
            param.coerce(1000)
        self.assertEqual(param.coerce("12"), 12)

    def test_messages_and_exports_carry_no_server_paths(self):
        folder = settings.UPLOAD_DIR / ("b" * 32) / "image" / "slide.btf"
        self.assertEqual(public_message(OSError(f"Unable to open {folder}")), "Unable to open image/slide.btf")
        uns = dict(deconvolution=dict(nnls=dict(log_path="/srv/x/run.log", cache_namespace="/srv/y", method="nnls")))
        self.assertEqual(_without_server_paths(uns), dict(deconvolution=dict(nnls=dict(method="nnls"))))


def _pool_by_loop(pooler):
    """The original implementation: one sparse row per polygon."""
    groups = pooler.filtered_adata.obs.groupby(['gdf_index'] + pooler.keep_cols, observed=True)
    counts = pooler.filtered_adata.X
    summed = sp.lil_matrix((groups.ngroups, counts.shape[1]))
    for row, (_, idx) in enumerate(groups.indices.items()):
        summed[row] = counts[idx].sum(0)
    return summed.tocsr()


FRAME_FIXTURE = Path(__file__).resolve().parents[1] / "frontend" / "src" / "api" / "__fixtures__" / "frame.bin"


def fixture_frame():
    """The frame decoded by frontend/src/api/frames.test.ts."""
    from hd_viewer.api import frames
    header = dict(layer="bins", base=0, version=1, row_version=3, rows=3, level="agg:2", constants=dict(width=2.0),
                  palette=["#ffffff", "#123456"], tooltips=[["Barcode", "$index{custom}"]], visible=True, filtered=True)
    ops = [("append", dict(x=np.array([1.5, 2.5, 3.5], dtype=np.float32), color=np.array([1, 0, 1], dtype=np.uint8),
                           codes=np.array([0, 300, 65535], dtype=np.uint16),
                           x_rings=[np.array([0.0, 1.0]), np.array([2.0]), np.array([3.0, 4.0, 5.0])])),
           ("set", "alpha", np.array([0.25, 0.5, 1.0], dtype=np.float32))]
    return frames.encode(header, ops, extra=dict(visible=np.array([1, 0, 1], dtype=np.uint8)))


class FrameFixtureTests(unittest.TestCase):
    def test_front_end_fixture_matches_the_encoder(self):
        """If the frame format changes, rewrite the fixture with
        `VHD_WRITE_FIXTURES=1 python -m unittest tests.test_server_hardening`."""
        import os
        frame = fixture_frame()
        if os.environ.get("VHD_WRITE_FIXTURES"):
            FRAME_FIXTURE.write_bytes(frame)
        self.assertEqual(FRAME_FIXTURE.read_bytes(), frame)


class ImagesDirTests(unittest.TestCase):
    def test_environment_variable_overrides_the_configured_image_folder(self):
        from hd_viewer.plugins.runners import ApptainerRunner, make_runner
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder) / "cellpose.sif").touch()
            plugin = SimpleNamespace(name="cellpose", image="ghcr.io/x/cellpose:1")
            config = dict(runner="apptainer", apptainer=dict(images_dir="images"))
            kwargs = dict(sdk_dir=Path(folder), segmentation_root=Path(folder) / "api")
            with patch.dict("os.environ", {"VHD_DECONV_IMAGES_DIR": folder}):
                seg = make_runner(config, **kwargs, runner_types={"apptainer": ApptainerRunner})
                deconv = make_runner(config, **kwargs, env_var="VHD_DECONV_RUNNER", kind="deconvolution",
                                     runner_types={"apptainer": ApptainerRunner})
            self.assertEqual(seg.image_uri(plugin), "docker://ghcr.io/x/cellpose:1")   # VHD_SEG_IMAGES_DIR unset
            self.assertEqual(deconv.image_uri(plugin), str(Path(folder) / "cellpose.sif"))
            self.assertEqual(config["apptainer"]["images_dir"], "images")               # config not modified


class InitialAggregationTests(unittest.TestCase):
    def test_empty_margins_do_not_coarsen_the_bins(self):
        """An 8 µm slide: 838 bins across inside an image three times wider."""
        from hd_viewer.core import viewport
        extent = (0.0, 838.0, 0.0, 838.0)                       # in native bins
        window = (-1500.0, 2400.0, -500.0, 1400.0)               # whole image, screen-shaped
        width, height = viewport.window_size(*viewport.clip_to_extent(*window, extent))
        choose = lambda span: viewport.aggregation_factor(span, (1, 2, 4, 8, 12), 1, target_bins_per_axis=360,
                                                          hysteresis_ratio=0.15)[0]
        self.assertEqual(choose(max(width, height)), 4)
        self.assertEqual(choose(max(viewport.window_size(*window))), 12)  # the margins counted
        self.assertEqual(viewport.clip_to_extent(*window, (5000, 6000, 0, 1)), window)   # no overlap
        self.assertEqual(viewport.clip_to_extent(*window, None), window)

    def test_bins_extent_covers_whole_edge_bins(self):
        import geopandas as gpd
        from shapely.geometry import box
        from hd_viewer.startup.sample import Bins
        gdf = gpd.GeoDataFrame({"index": ["a", "b"]}, geometry=[box(0, 0, 2, 2), box(8, 4, 10, 6)])
        adata = ad.AnnData(np.ones((2, 1)), obs=pd.DataFrame(index=["a", "b"]))
        sample = Bins(gdf, adata, pixels_per_bin=2, bin_size=8)
        self.assertEqual(sample.data_extent(), (-1.0, 11.0, -1.0, 7.0))


class RCTDCoordinateTests(unittest.TestCase):
    def test_coordinates_are_the_centres_of_each_observation(self):
        import geopandas as gpd
        from shapely.geometry import box
        from hd_viewer.core.deconvolution import observation_coordinates, snapshot_target
        import threading
        gdf = gpd.GeoDataFrame(dict(index=["b", "a", "a"]),
                               geometry=[box(0, 0, 2, 2), box(10, 10, 12, 12), box(20, 10, 22, 12)])
        adata = ad.AnnData(X=sp.csr_matrix(np.ones((2, 3), dtype=np.float32)), obs=pd.DataFrame(index=["a", "b"]))
        sample = SimpleNamespace(adata=adata, gdf=gdf, on="index", fetch_lock=threading.Lock(),
                                 _compute_centroids=lambda g: (g.geometry.centroid.x.to_numpy(), g.geometry.centroid.y.to_numpy()))
        np.testing.assert_allclose(observation_coordinates(sample), [[16, 11], [1, 1]])
        np.testing.assert_allclose(snapshot_target(sample, "spacexr").obsm["spatial"], [[16, 11], [1, 1]])
        self.assertNotIn("spatial", snapshot_target(sample, "nnls").obsm)
        sample.gdf = gdf.iloc[:1]                                   # "a" has no geometry
        self.assertIsNone(observation_coordinates(sample))


class FastPathTests(unittest.TestCase):
    def test_pooling_by_matrix_product_matches_the_loop(self):
        rng = np.random.default_rng(0)
        n, genes = 400, 30
        X = sp.random(n, genes, density=0.2, format="csr", random_state=1, data_rvs=lambda k: rng.integers(1, 9, k))
        obs = pd.DataFrame(dict(gdf_index=rng.integers(0, 37, n), array_row=0, array_col=0),
                           index=[f"b{i}" for i in range(n)])
        obs.loc[obs.index[:20], "gdf_index"] = np.nan                # bins in no polygon
        for keep in ([], ["array_row", "array_col"]):
            pooler = Pooler.__new__(Pooler)
            pooler.keep_cols = keep
            pooler.filtered_adata = ad.AnnData(X=X, obs=obs.copy(), var=pd.DataFrame(index=[f"g{j}" for j in range(genes)]))
            pooled = pooler.pool()
            np.testing.assert_allclose(pooled.X.toarray(), _pool_by_loop(pooler).toarray())
            self.assertEqual(pooled.n_obs, pooler.filtered_adata.obs["gdf_index"].nunique())

    def test_batched_de_matches_one_test_per_set(self):
        rng = np.random.default_rng(0)
        X = sp.random(300, 40, density=0.3, format="csc", random_state=2, data_rvs=lambda k: rng.poisson(3, k) + 1)
        adata = ad.AnnData(X=X.astype(np.float32))
        sc.pp.log1p(adata)
        labels = rng.integers(0, 4, 300)
        masks = {name: labels == i for i, name in enumerate(["A", "B", "C"])}   # label 3: in no set
        batched = DEAnalysis._de_vs_rest_batched(adata, list(masks), masks)
        for name, mask in masks.items():
            adata.obs["_one"] = np.where(mask, name, "Rest")
            sc.tl.rank_genes_groups(adata, groupby="_one", groups=[name], reference="Rest", method=settings.DE_METHOD)
            single = sc.get.rank_genes_groups_df(adata, group=name).set_index("names").sort_index()
            together = batched[name].set_index("names").sort_index()
            np.testing.assert_allclose(together["scores"], single["scores"], rtol=1e-6)
            np.testing.assert_allclose(together["pvals_adj"], single["pvals_adj"], rtol=1e-6)
            np.testing.assert_allclose(together["logfoldchanges"], single["logfoldchanges"], rtol=1e-5)

    def test_overlapping_sets_are_not_batched(self):
        adata = ad.AnnData(X=sp.csc_matrix(np.ones((10, 3), dtype=np.float32)))
        masks = {"A": np.arange(10) < 6, "B": np.arange(10) > 3}
        self.assertIsNone(DEAnalysis._de_vs_rest_batched(adata, ["A", "B"], masks))


if __name__ == "__main__":
    unittest.main()
