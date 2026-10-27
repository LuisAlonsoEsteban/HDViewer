"""Native Bokeh views of the plot data prepared by DEPresenter."""

from bokeh.layouts import column
from bokeh.models import ColumnDataSource, Div, FixedTicker, HoverTool, Legend, LegendItem
from bokeh.plotting import figure


def build_de_plots(plot_data):
    """Build fresh models on the document thread; no image files or HTML embeds."""
    plots = []
    for data in plot_data:
        p = figure(title=f"{data.name} vs {data.rest_label} · top genes", height=340,
                   sizing_mode="stretch_width", min_width=320,
                   tools="pan,wheel_zoom,box_zoom,box_select,tap,reset,save",
                   active_scroll="wheel_zoom", background_fill_color="#0a0a0a",
                   border_fill_color="#0a0a0a", outline_line_color="#444444")
        p.xaxis.ticker = FixedTicker(ticks=list(range(len(data.genes))))
        p.xaxis.major_label_overrides = dict(enumerate(data.genes))
        p.xaxis.major_label_orientation = .4
        p.yaxis.axis_label = "Log2 fold change" if data.kind == "dot" else "Expression · log1p(counts)"
        p.title.text_color = "#dddddd"
        p.axis.major_label_text_color = "#cccccc"
        p.axis.axis_label_text_color = "#cccccc"
        p.axis.axis_line_color = "#555555"
        p.grid.grid_line_color = "#333333"
        p.grid.grid_line_alpha = .4
        p.toolbar.logo = None
        hover_renderers, legend_items = [], []
        for group in dict.fromkeys(row["group"] for row in data.rows):
            rows = [row for row in data.rows if row["group"] == group]
            source = ColumnDataSource({key: [row[key] for row in rows] for key in rows[0]})
            if data.kind == "dot":
                dots = p.scatter("x", "lfc", source=source, size="size", color="color",
                                 alpha=.85, line_color="#eeeeee", selection_line_color="#ffffff",
                                 nonselection_alpha=.2, name="de_dots")
                hover_renderers.append(dots)
                legend_items.append(LegendItem(label=group, renderers=[dots]))
            else:
                violin = p.patches("xs", "ys", source=source, color="color", fill_alpha=.65,
                                   line_alpha=.9, nonselection_fill_alpha=.15, name="de_violins")
                quartiles = p.segment("x", "q1", "x", "q3", source=source,
                                      line_width=5, color="#dddddd")
                medians = p.scatter("x", "median", source=source, size=6, color="#ffffff")
                hover_renderers.extend([violin, medians])
                legend_items.append(LegendItem(label=group, renderers=[violin, quartiles, medians]))
        tooltips = [("Gene", "@gene"), ("Group", "@group"), ("Comparison", "@comparison"),
                    ("Log2 fold change", "@lfc{0.000}"), ("Adjusted p-value", "@padj"),
                    ("Expressing", "@pct{0.0}%"), ("Observations", "@n{0,0}"),
                    ("Mean log1p expression", "@mean{0.000}")]
        if data.kind == "violin":
            tooltips += [("Values in violin", "@n_shown{0,0}"), ("Q1", "@q1{0.000}"),
                         ("Median", "@median{0.000}"), ("Q3", "@q3{0.000}")]
        p.add_tools(HoverTool(renderers=hover_renderers, tooltips=tooltips))
        p.add_layout(Legend(items=legend_items, click_policy="hide", label_text_color="#cccccc",
                            background_fill_color="#111111", border_line_color="#444444"), "below")
        plots.append(p)
    return column(*plots, sizing_mode="stretch_width", spacing=16) if plots else column(
        Div(text="No positive marker genes to plot."), sizing_mode="stretch_width")
