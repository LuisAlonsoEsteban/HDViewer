import unittest

from bokeh.models import Div, Range1d
from bokeh.plotting import figure

from hd_viewer.board.layout.thumbnail import ViewportReadout


class ViewportReadoutTests(unittest.TestCase):
    def setUp(self):
        self.figure = figure(x_range=(0, 1000), y_range=(0, 1000))
        self.readout = ViewportReadout(self.figure)

    def test_plot_units_until_the_scale_is_known(self):
        self.assertEqual(self.readout.bbox(), (0, 0, 1000, 1000))

    def test_full_res_pixels(self):
        # Image 2000 px wide and 1000 px high stretched to a 1000-unit square
        self.readout.set_scale(2.0, 1.0, um_per_px=0.5)
        self.assertEqual(self.readout.bbox(), (0, 0, 2000, 1000))

    def test_follows_replaced_ranges(self):
        self.readout.set_scale(2.0, 1.0, um_per_px=0.5)
        self.figure.x_range = Range1d(100, 300)
        self.assertEqual(self.readout.bbox()[::2], (200, 600))

    def test_only_a_copy_button_with_no_live_updates(self):
        self.assertEqual(self.readout.copy_button.label, "Copy viewport coordinates")
        self.assertIsNotNone(self.readout.copy_button.icon)
        self.assertFalse(hasattr(self.readout, "text"))
        for axis in (self.figure.x_range, self.figure.y_range):
            self.assertFalse(axis.js_property_callbacks)

    def test_shown_in_the_top_right_of_another_plot(self):
        self.assertIsNone(self.readout.panel)
        thumbnail = figure(x_range=(0, 1000), y_range=(0, 1000))
        self.readout.place_on(thumbnail)
        self.assertIn(self.readout.panel, thumbnail.elements)
        self.assertNotIn(self.readout.panel, self.figure.elements)
        self.assertEqual(self.readout.panel.elements, [self.readout.copy_button])
        self.assertFalse(any(isinstance(e, Div) for e in self.readout.panel.elements))
        position = self.readout.panel.position
        self.assertEqual((position.target, position.symbol, self.readout.panel.anchor),
                         ("frame", "top_right", "top_right"))
        # It copies the main plot's viewport, not the one it is shown on
        self.figure.x_range.end = 500
        self.assertEqual(self.readout.bbox()[2], 500)


if __name__ == "__main__":
    unittest.main()
