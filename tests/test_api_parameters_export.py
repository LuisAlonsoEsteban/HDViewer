"""Parameters tab (extra glyph attributes), h5ad export and image renderer names."""
import io
import unittest

import anndata as ad
import numpy as np

from hd_viewer.api import frames
from tests.test_api import ApiTestCase, _objects


class ParametersAndExportTests(ApiTestCase):
    def setUp(self):
        super().setUp()
        self.client.get("/api/load")
        self.load_example()
        self.size = self.client.get("/api/workspace").json()["image"]["size"]
        self.client.post("/api/layers/bins/viewport",
                         json=dict(seq=1, x0=0, x1=self.size, y0=0, y1=self.size, width=600, height=600))
        session = next(iter(self.registry._sessions.values()))
        self.sample = session.views["bins"].sample
        obs = self.sample.adata.obs
        obs["region"] = np.where(np.arange(len(obs)) % 2, "a", "b").astype(object)

    def names(self, item):
        return {ref.get("attributes", {}).get("name") for ref in _objects(item)}

    def test_parameters_validate_change_the_figure_and_reach_the_frames(self):
        listed = self.client.get("/api/layers/bins/parameters").json()
        self.assertEqual([entry["attr"] for entry in listed], ["fill_alpha", "line_width", "line_color", "line_alpha"])
        self.assertTrue(all(entry["current"] is None for entry in listed))

        colours = self.client.get("/api/layers/bins/parameters/line_color/columns").json()
        self.assertIn("region", colours)
        self.assertNotIn(self.sample.adata.var_names[0], colours)
        widths = self.client.get("/api/layers/bins/parameters/line_width/columns").json()
        self.assertIn(self.sample.adata.var_names[0], widths)
        self.assertNotIn("region", widths)

        bad = [("line_width", dict(mode="static", value=99)),                    # out of range
               ("line_width", dict(mode="data", column="region")),               # not numeric
               ("line_color", dict(mode="data", column=self.sample.adata.var_names[0])),
               ("color", dict(mode="static", value="#fff")),                     # a core attribute
               ("line_alpha", dict(mode="other"))]
        for attr, body in bad:
            self.assertEqual(self.client.put(f"/api/layers/bins/parameters/{attr}", json=body).status_code, 422,
                             (attr, body))

        width = self.client.put("/api/layers/bins/parameters/line_width", json=dict(mode="static", value=3)).json()
        self.assertEqual((width["mode"], width["value"]), ("static", 3.0))
        gene = self.sample.adata.var_names[0]
        alpha = self.client.put("/api/layers/bins/parameters/line_alpha",
                                json=dict(mode="data", column=gene, min=0.2, max=0.9, hovered=True)).json()
        self.assertEqual((alpha["column"], alpha["min"], alpha["max"]), (gene, 0.2, 0.9))
        line = self.client.put("/api/layers/bins/parameters/line_color",
                               json=dict(mode="data", column="region", hovered=True)).json()
        self.assertEqual(line["column"], "region")

        # The rebuilt figure binds the glyph to the new columns, keeping the view
        item = self.client.get("/api/figures/map", params=dict(x0=10, x1=20, y0=5, y1=15)).json()
        self.assertIn("bins:mapper:line_color", self.names(item))
        text = str(item)
        self.assertIn("'field': 'line_width'", text)
        self.assertIn("'field': 'line_alpha'", text)

        header, ops, _ = frames.decode(self.client.get("/api/layers/bins/state").content)
        columns = ops[0][1]
        for name in ("line_width", "line_alpha", "line_color"):
            self.assertIn(name, columns)
        sequenced = np.asarray(columns["line_alpha"])
        self.assertTrue(((sequenced >= 0.2 - 1e-6) | np.isclose(sequenced, 0.5)).all())
        self.assertEqual(len(header["palettes"]["line_color"]), 3)              # not sequenced + a, b
        hovered = [row[0] for row in header["tooltips"]]
        self.assertIn("region", hovered)

        self.assertEqual(self.client.delete("/api/layers/bins/parameters/line_color").status_code, 204)
        header, ops, _ = frames.decode(self.client.get("/api/layers/bins/state").content)
        self.assertNotIn("line_color", ops[0][1])
        self.assertNotIn("line_color", header["palettes"])
        self.assertNotIn("bins:mapper:line_color", self.names(self.client.get("/api/figures/map").json()))
        self.assertEqual(self.client.delete("/api/layers/bins/parameters/nope").status_code, 404)

    def test_image_renderer_is_named(self):
        names = self.names(self.client.get("/api/figures/map").json())
        self.assertTrue({"image:url", "image:tiles"} & names)

    def test_export_streams_the_layer_with_its_results_and_removes_the_file(self):
        response = self.client.get("/api/layers/bins/export")
        self.assertEqual(response.status_code, 200, response.text)
        self.assertIn('filename="bins.h5ad"', response.headers["content-disposition"])
        exported = ad.read_h5ad(io.BytesIO(response.content))
        self.assertEqual(exported.shape, self.sample.adata.shape)
        self.assertEqual(list(exported.obs["region"].astype(str)[:2]), ["b", "a"])
        session = next(iter(self.registry._sessions.values()))
        self.assertEqual(list((session.directory / "exports").iterdir()), [])
        self.assertEqual(self.client.get("/api/layers/nope/export").status_code, 404)


if __name__ == "__main__":
    unittest.main()


class NonAggregatedSelectionTests(unittest.TestCase):
    def test_spots_and_cells_selections_map_to_rows(self):
        """Layers without aggregation (spots, cells) have no aggregate name prefix."""
        from types import SimpleNamespace
        import pandas as pd
        from hd_viewer.board.components.layer.utils.selection import native_selection, selection_indices
        spots = SimpleNamespace(adata=SimpleNamespace(obs_names=pd.Index(["s0", "s1", "s2"])))
        self.assertEqual(native_selection(spots, ["s1", "x"]), {"s1"})
        self.assertEqual(selection_indices(spots, {"s2", "s0"}, ["s0", "s1", "s2"]), [0, 2])


class SpotsKeepTheirRadiusTests(unittest.TestCase):
    def test_radius_survives_clearing_and_recolouring(self):
        """The circle glyph reads `radius` per row: no frame may drop it."""
        from hd_viewer.api.views import LayerView
        from tests.test_roi_manager import make_layer
        sample = make_layer("spots").sample
        view = LayerView("spots", "spots", sample, "Spots", "#ffffff")
        view.update_viewport(-1, 5, -1, 5, 600, 600)
        self.assertIn("radius", view.rows.sink.data)
        sample.adata.obs["region"] = ["a", "b"] * 8
        for ops in (view.colour_by_column("region"), view.colour_by_column("g1"),
                    view.set_parameter("line_width", mode="static", value=2) and view.rows.sink.take()):
            data = view.rows.sink.data
            self.assertIn("radius", data)
            lengths = {name: len(values) for name, values in data.items()}
            self.assertEqual(len(set(lengths.values())), 1, lengths)
        view.rows.clear_sources()
        self.assertIn("radius", view.rows.sink.data)
