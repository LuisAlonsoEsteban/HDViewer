import numpy as np
import scanpy as sc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Div, TabPanel, Tabs, DataTable, TableColumn, NumberFormatter
from bokeh.plotting import figure
from . import Panel
import math


class MetadataPanel(Panel):

    def __init__(self):
        super().__init__(name="metadata")
        # Mantener este diccionario para compatibilidad interna de la app
        self._marker_table_sources = {}
        self._marker_divs = {}
        self._de_divs = {}

    def _compute_general_stats(self, geom):
        """Calcula las estadísticas numéricas para la tabla de resumen."""
        adata = geom.sample.adata
        n_total = geom.sample.n_geometries
        n_not_seq = geom.sample.n_not_sequenced
        n_seq = n_total - n_not_seq
        pct_seq = round(100 * n_seq / n_total, 1) if n_total > 0 else 0

        mean_umis = 0.0
        if adata is not None and "total_counts" in adata.obs.columns:
            mean_umis = round(float(adata.obs["total_counts"].mean()), 1)
        
        return {
            "Values": [str(n_total), str(n_seq), str(n_not_seq), f"{pct_seq}%", str(mean_umis)]
        }

    def _build_general_tab(self, geometries):
        """Construye la pestaña 'General Info' con la tabla resumen en HTML."""
        rows_html = ""
        for geom in geometries:
            stats = self._compute_general_stats(geom)
            vals = stats["Values"]
            rows_html += f"""
            <tr>
                <td style="padding: 6px 10px; border-bottom: 1px solid #333;"><b>{geom.alias}</b></td>
                <td style="padding: 6px 10px; border-bottom: 1px solid #333;">{vals[0]}</td>
                <td style="padding: 6px 10px; border-bottom: 1px solid #333;">{vals[1]}</td>
                <td style="padding: 6px 10px; border-bottom: 1px solid #333;">{vals[2]}</td>
                <td style="padding: 6px 10px; border-bottom: 1px solid #333;">{vals[3]}</td>
                <td style="padding: 6px 10px; border-bottom: 1px solid #333;">{vals[4]}</td>
            </tr>"""

        html_content = f"""
        <div style="padding: 10px; color: #dddddd; font-family: sans-serif;">
            <h3 style="color: #dddddd; margin: 0 0 10px 0;">Dataset Summary</h3>
            <table style="width: 100%; border-collapse: collapse; font-size: 12px; text-align: left;">
                <thead>
                    <tr style="background: #2a2a3e; border-bottom: 1px solid #444;">
                        <th style="padding: 6px 10px; color: #5cb85c;">Sample</th>
                        <th style="padding: 6px 10px; color: #5cb85c;">Total</th>
                        <th style="padding: 6px 10px; color: #5cb85c;">Sequenced</th>
                        <th style="padding: 6px 10px; color: #5cb85c;">Not Seq</th>
                        <th style="padding: 6px 10px; color: #5cb85c;">% Seq</th>
                        <th style="padding: 6px 10px; color: #5cb85c;">Mean UMIs</th>
                    </tr>
                </thead>
                <tbody>{rows_html}</tbody>
            </table>
        </div>
        """
        return TabPanel(
            child=column(Div(text=html_content), sizing_mode="stretch_width"),
            title="General Info"
        )

    def populate(self, geometries):
        print("--- CONFIGURANDO INTERFAZ DE METADATOS ---")
        
        for geom in geometries:
            geom._metadata_panel = self

        # Tab principal: General Info
        #main_tabs = [self._build_general_tab(geometries)]
        main_tabs = []
        # Un tab principal por cada geom con subtabs QC y Markers dentro
        for geom in geometries:
            adata = geom.sample.adata
            if adata is None:
                continue

            try:
                adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
                if "pct_counts_mt" not in adata.obs.columns:
                    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=False)
            except Exception as e:
                print(f"Error QC para {geom.alias}: {e}")
                continue

            obs = adata.obs
            vals_mt = np.array(obs["pct_counts_mt"].dropna(), dtype=float)
            vals_genes = np.array(obs["n_genes_by_counts"].dropna(), dtype=float)

            # --- FIGURA 1: Histograma % MT ---
            hist, edges = np.histogram(vals_mt, bins=30)
            hist_src = ColumnDataSource(dict(
                top=hist.tolist(),
                left=edges[:-1].tolist(),
                right=edges[1:].tolist()
            ))
            p_hist_mt = figure(
                title="% MT", width=210, height=190,
                x_axis_label="% MT", y_axis_label="N",
                tools="reset,pan,wheel_zoom"
            )
            p_hist_mt.xaxis.major_label_orientation = 0.8  # rotar etiquetas
            p_hist_mt.xaxis.ticker.desired_num_ticks = 5   # menos ticks

            p_hist_mt.quad(
                top="top", bottom=0, left="left", right="right", source=hist_src,
                fill_color="#f14e4e", line_color="#1e1e2e", alpha=0.8
            )

            # --- FIGURA 2: Histograma Genes per spot ---
            hist2, edges2 = np.histogram(vals_genes, bins=30)
            hist_src2 = ColumnDataSource(dict(
                top=hist2.tolist(),
                left=edges2[:-1].tolist(),
                right=edges2[1:].tolist()
            ))
            p_hist_genes = figure(
                title="Genes per spot", width=210, height=190,
                x_axis_label="N genes", y_axis_label="Frecuencia",
                tools="reset,pan,wheel_zoom"
            )
            p_hist_genes.xaxis.major_label_orientation = 0.8  # rotar etiquetas
            p_hist_genes.xaxis.ticker.desired_num_ticks = 5   # menos ticks

            p_hist_genes.quad(
                top="top", bottom=0, left="left", right="right", source=hist_src2,
                fill_color="#4e9af1", line_color="#1e1e2e", alpha=0.8
            )

            # --- FIGURA 3: Scatter UMIs vs % MT ---
            scatter_src = ColumnDataSource(dict(
                total_counts=obs["total_counts"].tolist(),
                pct_counts_mt=obs["pct_counts_mt"].tolist()
            ))
            p_scatter_mt = figure(
                title="UMIs vs % MT", width=210, height=190,
                x_axis_label="UMIs", y_axis_label="% MT",
                tools="reset,pan,wheel_zoom"
            )
            p_scatter_mt.xaxis.major_label_orientation = 0.8  # rotar etiquetas
            p_scatter_mt.xaxis.ticker.desired_num_ticks = 5   # menos ticks

            p_scatter_mt.circle(
                "total_counts", "pct_counts_mt", source=scatter_src,
                size=3, color="#4ef1a1", alpha=0.6
            )

            # --- FIGURA 4: Scatter UMIs vs Genes ---
            scatter_src2 = ColumnDataSource(dict(
                total_counts=obs["total_counts"].tolist(),
                n_genes_by_counts=obs["n_genes_by_counts"].tolist()
            ))
            p_scatter_genes = figure(
                title="UMIs vs Genes", width=210, height=190,
                x_axis_label="UMIs", y_axis_label="N genes",
                tools="reset,pan,wheel_zoom"
            )

            p_scatter_genes.xaxis.major_label_orientation = 0.8  # rotar etiquetas
            p_scatter_genes.xaxis.ticker.desired_num_ticks = 5   # menos ticks

            p_scatter_genes.circle(
                "total_counts", "n_genes_by_counts", source=scatter_src2,
                size=3, color="#f1a14e", alpha=0.6
            )

            qc_tab = TabPanel(
                child=column(
                    row(p_hist_mt, p_hist_genes, spacing=5),
                    row(p_scatter_mt, p_scatter_genes, spacing=5),
                    sizing_mode="stretch_width",
                    styles={"overflow-y": "auto", "max-height": "420px"}
                ),
                title="QC"
            )
            marker_tab = self._build_marker_tab(geom)

            # Subtabs dentro del tab del geom
            inner_tabs = Tabs(tabs=[qc_tab, marker_tab, self._build_de_tab(geom)])

            main_tabs.append(TabPanel(
                child=column(inner_tabs, sizing_mode="stretch_width"),
                title=geom.alias
            ))
        main_tabs.insert(0, self._build_general_tab(geometries))

        if main_tabs:
            self.set_children([Tabs(tabs=main_tabs)])
        else:
            self.set_children([Div(text="<p style='color:gray'>Esperando datos...</p>")])

        print("--- INTERFAZ ACTUALIZADA ---")

    def _build_marker_tab(self, geom):
        marker_div = Div(
            text="<p style='color:#888; font-size:11px; padding:10px;'>"
                "Run clustering to populate.</p>"
        )
        self._marker_divs[geom.alias] = marker_div  # ← guardar referencia
        self._marker_table_sources[geom.alias] = ColumnDataSource(
            dict(cluster=[], gene=[], score=[])
        )
        return TabPanel(
            child=column(marker_div),
            title="Markers"
        )

    def update_marker_table(self, geom):
        print(f"--- [ACTUALIZANDO TABLA] Extrayendo marcadores para {geom.alias} ---")
        markers = geom.sample.get_marker_genes(n_genes=10)
        if not markers or geom.alias not in self._marker_divs:
            print(f"Aviso: No se encontraron marcadores para {geom.alias}.")
            return
        rows_html = ""
        for cluster, gene_list in markers.items():
            for gene, score, pval, logfc in gene_list:
                # Formatear pval en notación científica
                pval_str = f"{pval:.2e}" if (pval and pval > 0 and not math.isnan(pval)) else "0.00e+00"
                # Colorear logFC en rojo si es negativo, verde si positivo
                logfc_color = "#4ef1a1" if logfc >= 0 else "#f14e4e"
                rows_html += f"""
                <tr>
                    <td style='padding:3px 8px;color:#aaa'>{cluster}</td>
                    <td style='padding:3px 8px;color:#fff'>{gene}</td>
                    <td style='padding:3px 8px;color:#4e9af1'>{score}</td>
                    <td style='padding:3px 8px;color:#aaa'>{pval_str}</td>
                    <td style='padding:3px 8px;color:{logfc_color}'>{logfc}</td>
                </tr>"""
        html = f"""
        <div style='padding:10px;overflow-y:auto;max-height:350px'>
        <h4 style='color:#dddddd;margin:0 0 8px 0;position:sticky;top:0;background:#1e1e2e'>
            Top Marker Genes</h4>
        <table style='border-collapse:collapse;font-size:12px;width:100%'>
            <thead style='position:sticky;top:28px;background:#1e1e2e'>
                <tr style='border-bottom:1px solid #444'>
                    <th style='padding:3px 8px;color:#888;text-align:left'>Cluster</th>
                    <th style='padding:3px 8px;color:#888;text-align:left'>Gene</th>
                    <th style='padding:3px 8px;color:#888;text-align:left'>Z-Score</th>
                    <th style='padding:3px 8px;color:#888;text-align:left'>Adj P-val</th>
                    <th style='padding:3px 8px;color:#888;text-align:left'>Log2FC</th>
                </tr>
            </thead>
            <tbody>{rows_html}</tbody>
        </table></div>"""
        self._marker_divs[geom.alias].text = html
        print(f"--- [TABLA ACTUALIZADA] ---")
    
    def _build_de_tab(self, geom):
        de_div = Div(
            text="<p style='color:#888; font-size:11px; padding:10px;'>"
                "Create groups and run DE to populate.</p>"
        )
        self._de_divs[geom.alias] = de_div
        return TabPanel(
            child=column(de_div, sizing_mode="stretch_width"),
            title="Groups DE"
        )
    
    def update_de_table(self, geom, cname, result_html):
        print(f"[update_de_table] alias={geom.alias}, in_de_divs={geom.alias in self._de_divs}")
        if geom.alias in self._de_divs:
            print(f"[update_de_table] Actualizando div, html len={len(result_html)}")
            #self._de_divs[geom.alias].text = result_html
            #result_html = """
            #<div style="
            #background:red;
            #height:400px;
            #width:100%;
            #font-size:50px;
            #color:white">
            #PLOT TEST
            #</div>
            #"""
            self._de_divs[geom.alias].text = """
            <h1 style="color:red">HELLO</h1>
            <div style="background:red;height:300px;width:300px">
            TEST
            </div>
            """

            print(type(self._de_divs[geom.alias]))
            print(f"[update_de_table] Div actualizado")