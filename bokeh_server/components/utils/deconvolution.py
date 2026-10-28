"""Per-layer reference loading, background analysis and spatial result display."""

from ...backend_config import deconvolution_kwargs

from ... import settings

from html import escape

from bokeh.models import Button, Div, Select

from src.utils.deconvolution import deconvolve
from .loading import LoadingLogger


class DeconvolutionManager:

    DECONV_RESULT_PREFIX = settings.DECONV_RESULT_PREFIX

    def __init__(self, component):
        self.component = component
        self.reference = None
        self.reference_path = None
        self.reference_name = None
        self.label_column = ""
        self.columns = []

    def apply_result(self, result, label_column):
        fractions, genes = result
        sample = self.component.sample
        if not fractions.index.equals(sample.adata.obs_names):
            raise ValueError("Layer observations changed during deconvolution. Run it again.")

        run = 1
        while any(str(c).startswith(f"{self.DECONV_RESULT_PREFIX}{run}:") for c in sample.get_colnames()):
            run += 1
        columns = [f"{self.DECONV_RESULT_PREFIX}{run}:{name}" for name in fractions.columns]
        # The dominant-type label makes the whole NumPy matrix object-typed.
        # Preserve each column's dtype so fractions stay continuous on the map.
        numeric = set(fractions.select_dtypes(include="number").columns)
        fraction_columns = [column for name, column in zip(fractions.columns, columns) if name in numeric]
        with sample.fetch_lock:
            for original, name in zip(fractions.columns, columns):
                values = fractions[original]
                sample.adata.obs[name] = values.astype(float) if original in numeric else values.astype("category")
            sample.adata.obsm["deconvolution"] = fractions
            sample.adata.uns["deconvolution"] = {
                "method": "nnls", "reference_path": self.reference_path or "",
                "reference_name": self.reference_name or "",
                "label_column": label_column, "cell_types": [name for name in fractions.columns if name in numeric],
                "genes": genes, "fraction_columns": fraction_columns,
            }
            # Numeric metadata normally sums at coarse bin levels; fractions average.
            means = list(sample.adata.uns.get("mean_obs_columns", []))
            sample.adata.uns["mean_obs_columns"] = means + fraction_columns
            if hasattr(sample, "clear_multiscale_cache"):
                sample.clear_multiscale_cache()
        self.columns = fraction_columns
        if self.component.autocomplete_input is not None:
            self.component.autocomplete_input.completions = self.component.get_colnames()
        if hasattr(self.component, "add_param_form"):
            self.component.add_param_form.refresh_completions()

    def get_widgets(self, doc):
        options = [] if self.reference is None else [str(c) for c in self.reference.obs.columns]
        labels = Select(title="Reference cell-type column (.obs)", options=options, disabled=not options,
                        value=self.label_column if self.label_column in options else (options[0] if options else ""))
        labels.on_change("value", lambda attr, old, new: setattr(self, "label_column", new))
        run = Button(label="Run deconvolution", button_type="success", disabled=not options)
        display = Select(title="Display cell-type fraction", options=self.columns,
                         value=self.columns[0] if self.columns else "", disabled=not self.columns)
        status = LoadingLogger(doc)
        if self.reference is None:
            status.set_message("No single-cell reference was provided for this dataset. "
                               "Deconvolution is unavailable; the other layer tools remain available.")
        elif not options:
            status.set_message("The reference contains counts but no .obs annotation columns. "
                               "Upload a barcode-indexed annotation CSV or Parquet file alongside it in the data loading step.")
        else:
            status.set_message(escape(
                f"Reference ready: {self.reference_name or 'loaded reference'} "
                f"({self.reference.n_obs:,} cells, {self.reference.n_vars:,} genes). "
                "Choose the cell-type column and run deconvolution."
            ))

        def on_run():
            if status.running or self.reference is None or not labels.value:
                return
            reference, label_column = self.reference, labels.value
            run.disabled = labels.disabled = True

            def restore():
                run.disabled = labels.disabled = False

            def done(result):
                try:
                    self.apply_result(result, label_column)
                    display.options = self.columns
                    display.disabled = False
                    display.value = self.columns[0]
                    show()
                    status.success(f"Deconvolved {len(result[0]):,} observations into {len(self.columns)} cell types")
                except Exception as exc:
                    status.error(exc)
                finally:
                    restore()

            def failed(exc):
                restore()
                status.error(exc)

            status.run("Deconvolution (NNLS)",
                       lambda: deconvolve(self.component.sample.adata, reference, label_column, **deconvolution_kwargs()), done, failed)

        def show():
            if not display.value:
                return
            if self.component.autocomplete_input is not None:
                self.component.autocomplete_input.value = display.value
            self.component.show_colname(display.value)
            if hasattr(self.component, "loader"):
                self.component.loader.force_reload()

        display.on_change("value", lambda attr, old, new: show())
        run.on_click(on_run)
        return [Div(text="NNLS estimates expression mixture fractions from an annotated single-cell "
                         "reference. Use linear counts in .X and matching gene identifiers. "
                         "Runs on all native observations of this layer; empty observations remain unassigned."),
                labels, run, status.div, display]
