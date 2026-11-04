from types import SimpleNamespace
import json
import unittest

import anndata as ad
from bokeh.document import Document
from bokeh.embed import json_item
from bokeh.models import BoxSelectTool, HoverTool, Legend, PanTool, Plot, ResetTool, WheelZoomTool
import numpy as np
import pandas as pd

from hd_viewer.board.operations.layer.de import DEResult
from hd_viewer.board.operations.layer.de import DEPresenter
from tests.test_de import log1p_normalized, make_analysis
from hd_viewer.board.operations.layer.de import build_de_plots


class DEPresentationTests(unittest.TestCase):
    def test_presentation_builds_native_interactive_figures(self):
        for alias, kind in (("Spot", "dot"), ("Cell", "violin"), ("Bin", "violin")):
            with self.subTest(alias=alias):
                analysis, masks = make_analysis()
                original = analysis.sample.adata.copy()
                result = analysis.run(masks, mode=1)
                view = DEPresenter(alias=alias, colors={"A": "#123456"}).prepare(result)
                self.assertIsInstance(result, DEResult)
                self.assertEqual(view.header, "DE: A vs B")
                top = {name: frame.sort_values("scores", ascending=False).names.iloc[0]
                       for name, frame in result.results.items()}
                self.assertEqual(top, {"A": "g0", "B": "g1"})
                self.assertEqual(len(view.plot_data), 2)
                self.assertTrue(all(data.kind == kind for data in view.plot_data))
                row = next(row for row in view.plot_data[0].rows if row['gene'] == 'g0' and row['group'] == 'A')
                expected = log1p_normalized(original.X)[masks['A'], 0]
                self.assertAlmostEqual(row['mean'], expected.mean())
                self.assertEqual(row['n'], 12)
                self.assertAlmostEqual(row['pct'], 100 * (expected > 0).mean())
                if kind == "violin":
                    self.assertAlmostEqual(row['median'], np.median(expected[expected > 0]))
                    self.assertTrue(np.isfinite(row['xs']).all())
                    self.assertEqual(len(row['xs']), len(row['ys']))
                layout = build_de_plots(view.plot_data)
                plots = list(layout.select({'type': Plot}))
                self.assertEqual(len(plots), 2)
                for plot in plots:
                    hover = next(tool for tool in plot.tools if isinstance(tool, HoverTool))
                    self.assertTrue(hover.renderers)
                    self.assertIsNotNone(plot.toolbar.active_inspect)
                    for tool_type in (PanTool, WheelZoomTool, BoxSelectTool, ResetTool):
                        self.assertTrue(any(isinstance(tool, tool_type) for tool in plot.tools))
                    self.assertTrue(all(legend.click_policy == 'hide' for legend in plot.select({'type': Legend})))
                    fields = dict(hover.tooltips)
                    self.assertEqual(fields['Gene'], '@gene')
                    self.assertEqual(fields['Group'], '@group')
                    self.assertIn('Adjusted p-value', fields)
                    for renderer in hover.renderers:
                        source = renderer.data_source.data
                        self.assertEqual(len({len(column) for column in source.values()}), 1)
                        self.assertTrue({'gene', 'group', 'lfc', 'padj', 'pct', 'n'} <= set(source))
                serialized = json.dumps(json_item(build_de_plots(view.plot_data)))
                self.assertIn('HoverTool', serialized)
                self.assertNotIn('data:image/', serialized)
                self.assertNotIn('<img', serialized)
                self.assertNotIn('static/data', serialized)
                document = Document()
                document.add_root(layout)
                document.validate()
                document.to_json()
                np.testing.assert_array_equal(analysis.sample.adata.X, original.X)
                pd.testing.assert_frame_equal(analysis.sample.adata.obs, original.obs)

    def test_constant_and_zero_distributions_have_finite_density_coordinates(self):
        analysis, masks = make_analysis()
        analysis.sample.adata.X[:] = 0
        analysis.sample.adata.X[:12, 0] = 20
        analysis.sample.adata.X[12:24, 1] = 20
        result = analysis.run(masks, mode=1)
        view = DEPresenter(alias="Cell").prepare(result)
        self.assertTrue(view.plot_data)
        for data in view.plot_data:
            for row in data.rows:
                self.assertTrue(np.isfinite(row['xs']).all())
                self.assertTrue(np.isfinite(row['ys']).all())
                self.assertEqual(len(row['xs']), len(row['ys']))
        Document().add_root(build_de_plots(view.plot_data))


if __name__ == "__main__":
    unittest.main()
