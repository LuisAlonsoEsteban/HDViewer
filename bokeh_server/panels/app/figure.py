from ... import settings

from . import Panel

from bokeh. plotting import figure

from bokeh.models import Div, CustomJS
from bokeh.layouts import column

from ...components.measure import MeasureTool
from ...components.roi import ROIReadout

class FigurePanel(Panel):

    p, img_alpha_slider, measure_tool, roi_readout = None, None, None, None

    def __init__(self):
        super().__init__(name="plot_column")
        self.layout.sizing_mode = "stretch_both"
        self.layout.min_height = 0
        self.de_layout = column(name="de_results", sizing_mode="stretch_width")

        self.image_attached = False
        self.alpha_setup = False
    
    def create_figure(self, max_dim):
        self.p = figure(
            sizing_mode='stretch_both',
            min_height=settings.PLOT_MIN_HEIGHT,
            x_range=(0, max_dim), 
            y_range=(0, max_dim),
            background_fill_color=settings.PLOT_BACKGROUND_COLOR,
            tools='pan,wheel_zoom,box_select,lasso_select,tap,save,reset')
        
        self.p.grid.visible = False
        self.p.xaxis.visible = False
        self.p.yaxis.visible = False
        self.p.toolbar.logo = None
        self.p.toolbar.autohide = True
        self.p.toolbar.active_inspect = None

        self._preserve_image_aspect()

        self.measure_tool = MeasureTool(self.p)
        self.roi_readout = ROIReadout(self.p)

    def _preserve_image_aspect(self):
        callback = CustomJS(args=dict(plot=self.p), code="""
            let scheduled = false
            let observer = null
            export default ({plot}) => {
                const schedule = () => {
                    if (scheduled) return
                    scheduled = true
                    requestAnimationFrame(() => {
                        try {
                            if (plot.document == null) { observer?.disconnect(); return }
                            const w = plot.inner_width, h = plot.inner_height
                            const xr = plot.x_range, yr = plot.y_range
                            const dx = xr.end-xr.start, dy = yr.end-yr.start
                            if (!(w > 0 && h > 0 && dx > 0 && dy > 0)) return
                            const units = Math.max(dx/w, dy/h)
                            const width = units*w, height = units*h
                            if (Math.abs(width-dx) <= dx*1e-7 && Math.abs(height-dy) <= dy*1e-7) return
                            const cx = (xr.start+xr.end)/2, cy = (yr.start+yr.end)/2
                            xr.setv({start:cx-width/2, end:cx+width/2})
                            yr.setv({start:cy-height/2, end:cy+height/2})
                        } finally { scheduled = false }
                    })
                }
                if (observer == null) {
                    const view = Bokeh.index.find_one(plot)
                    if (view != null) {
                        observer = new ResizeObserver(schedule)
                        observer.observe(view.el)
                    }
                }
                schedule()
            }
        """)
        self.p.js_on_event('document_ready', callback)
        for rng in (self.p.x_range, self.p.y_range):
            rng.js_on_change('start', callback)
            rng.js_on_change('end', callback)

    def attach_image(self, image_component, w=None, h=None):
        if self.p:
            image_component.attach_image(self.p, w=w, h=h)
            self.image_attached = True

    def setup_image_alpha(self, image_component):
        if self.p:
            image_component.add_img_overlay(self.p)
            self.img_alpha_slider = image_component.get_alpha_slider()
            self.alpha_setup = True
    
    def populate(self, image_component=None, w=None, h=None):
        if image_component is not None:
            if self.p is None:
                self.create_figure(image_component.image_processor.max_dim)
            self.set_children([self.p])

            if not self.image_attached:
                self.attach_image(image_component, w=w, h=h)
            if not self.alpha_setup:
                self.setup_image_alpha(image_component)

            self.add_children(self.img_alpha_slider)

    



        
