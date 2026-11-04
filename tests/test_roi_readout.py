import unittest

from bokeh.models import Range1d
from bokeh.plotting import figure

from hd_viewer.board.layout.thumbnail import ViewportReadout as ROIReadout


class ROIReadoutTests(unittest.TestCase):
    def setUp(self):
        self.figure = figure(x_range=(0, 1000), y_range=(0, 1000))
        self.readout = ROIReadout(self.figure)

    def test_plot_units_until_the_scale_is_known(self):
        self.assertEqual(self.readout.text.text, "x  0 – 1,000\ny  0 – 1,000")

    def test_full_res_pixels_and_physical_size(self):
        # Image 2000 px wide and 1000 px high stretched to a 1000-unit square
        self.readout.set_scale(2.0, 1.0, um_per_px=0.5)
        self.assertEqual(self.readout.bbox(), (0, 0, 2000, 1000))
        self.assertEqual(self.readout.text.text,
                         "x  0 – 2,000\ny  0 – 1,000\n1.00 mm × 500 µm")

    def test_follows_replaced_ranges(self):
        self.readout.set_scale(2.0, 1.0, um_per_px=0.5)
        self.figure.x_range = Range1d(100, 300)
        self.assertEqual(self.readout.bbox()[::2], (200, 600))
        self.assertIn("x  200 – 600", self.readout.text.text)
        # The browser-side update is attached to the new range too
        self.assertIn(self.readout._update, self.figure.x_range.js_property_callbacks["change:start"])

    def test_shown_in_the_top_right_of_another_plot(self):
        self.assertIsNone(self.readout.panel)
        thumbnail = figure(x_range=(0, 1000), y_range=(0, 1000))
        self.readout.place_on(thumbnail)
        self.assertIn(self.readout.panel, thumbnail.elements)
        self.assertNotIn(self.readout.panel, self.figure.elements)
        position = self.readout.panel.position
        self.assertEqual((position.target, position.symbol, self.readout.panel.anchor),
                         ("frame", "top_right", "top_right"))
        # It keeps tracking the main plot, not the one it is shown on
        self.figure.x_range.end = 500
        self.readout.render()
        self.assertIn("x  0 – 500", self.readout.text.text)
        self.assertIsNotNone(self.readout.copy_button.icon)
        self.assertEqual(self.readout.copy_button.label, "")


if __name__ == "__main__":
    unittest.main()
