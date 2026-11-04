from types import SimpleNamespace
import unittest

import anndata as ad
import numpy as np
import pandas as pd

from hd_viewer.board.operations.layer.de import DEAnalysis, DEResult


def log1p_normalized(X, target_sum=1e4):
    """What DE tests: log1p of counts scaled to `target_sum` per observation."""
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float)
    return np.log1p(X / X.sum(axis=1, keepdims=True) * target_sum)


def make_analysis():
    rng = np.random.default_rng(12)
    counts = rng.poisson(1., (36, 12)).astype(float)
    counts[:12, 0] += 30
    counts[12:24, 1] += 30
    data = ad.AnnData(counts, obs=pd.DataFrame(index=[f"s{i}" for i in range(36)]),
                      var=pd.DataFrame(index=[f"g{i}" for i in range(12)]))
    masks = {"A": np.arange(36) < 12, "B": (np.arange(36) >= 12) & (np.arange(36) < 24)}
    return DEAnalysis(SimpleNamespace(adata=data)), masks


class IndependentDETests(unittest.TestCase):
    def test_analysis_returns_statistics_and_aligned_data_without_presentation(self):
        analysis, masks = make_analysis()
        original = analysis.sample.adata.copy()
        result = analysis.run(masks, mode=1)
        self.assertIsInstance(result, DEResult)
        top = {name: frame.sort_values("scores", ascending=False).names.iloc[0]
               for name, frame in result.results.items()}
        self.assertEqual(top, {"A": "g0", "B": "g1"})
        self.assertEqual(result.mode, 1)
        self.assertEqual(result.names, ("A", "B"))
        self.assertEqual(result.expression.n_obs, 24)
        self.assertEqual(result.expression.obs_names.tolist(), original.obs_names[:24].tolist())
        np.testing.assert_allclose(result.expression.X, log1p_normalized(original.X)[:24])
        self.assertEqual(int(result.masks['A'].sum()), 12)
        self.assertEqual(int(result.masks['B'].sum()), 12)
        self.assertFalse((result.masks['A'] & result.masks['B']).any())
        for presentation_field in ('header', 'merged', 'table_html', 'plot_data', 'colors'):
            self.assertFalse(hasattr(result, presentation_field), presentation_field)
        np.testing.assert_array_equal(analysis.sample.adata.X, original.X)
        pd.testing.assert_frame_equal(analysis.sample.adata.obs, original.obs)

    def test_result_memberships_survive_later_runs_and_input_changes(self):
        analysis, masks = make_analysis()
        result = analysis.run(masks)
        expected = result.masks['A'].copy()
        masks['A'][:] = False
        analysis.run({'B': masks['B']})
        np.testing.assert_array_equal(result.masks['A'], expected)
        self.assertEqual(result.expression.obs.columns.tolist(), [])
        self.assertEqual(result.expression.n_obs, 36)

    def test_overlap_is_excluded_from_result_expression_and_memberships(self):
        analysis, masks = make_analysis()
        masks['A'][12:15] = True
        result = analysis.run(masks, mode=1)
        self.assertEqual(result.n_overlap, 3)
        self.assertEqual(result.expression.n_obs, 21)
        self.assertFalse(set(analysis.sample.adata.obs_names[12:15]) & set(result.expression.obs_names))
        self.assertEqual(int(result.masks['A'].sum()), 12)
        self.assertEqual(int(result.masks['B'].sum()), 9)

    def test_masks_and_modes_are_validated_before_analysis(self):
        analysis, masks = make_analysis()
        for invalid, mode in (({}, 0), ({"A": [True]}, 0), (masks, 2)):
            with self.subTest(mode=mode, names=list(invalid)):
                with self.assertRaises(ValueError):
                    analysis.run(invalid, mode)


if __name__ == "__main__":
    unittest.main()
