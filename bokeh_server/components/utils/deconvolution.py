"""Per-layer reference loading, background analysis and spatial result display."""

from pathlib import Path

import anndata as ad
from bokeh.models import Button, Div, Select, TextInput

from src.deconvolution import deconvolve
from .loading import LoadingLogger


class DeconvolutionManager:
    def __init__(self, component):
        self.component = component
        self.reference = None
        self.reference_path = None
        self.label_column = ""
        self.columns = []

    def apply_result(self, result, label_column):
        fractions, genes = result
        sample = self.component.sample
        if not fractions.index.equals(sample.adata.obs_names):
            raise ValueError("Layer observations changed during deconvolution. Run it again.")
        # Versioned names keep previous runs usable by drawing parameters/groups.
        run = 1
        while any(str(c).startswith(f"deconv_{run}:") for c in sample.get_colnames()):
            run += 1
        columns = [f"deconv_{run}:{name}" for name in fractions.columns]
        with sample.fetch_lock:
            for name, values in zip(columns, fractions.to_numpy().T):
                sample.adata.obs[name] = values
            sample.adata.obsm["deconvolution"] = fractions
            sample.adata.uns["deconvolution"] = {
                "method": "nnls", "reference_path": self.reference_path or "",
                "label_column": label_column, "cell_types": fractions.columns.tolist(),
                "genes": genes, "fraction_columns": columns,
            }
            # Numeric metadata normally sums at coarse bin levels; fractions average.
            means = list(sample.adata.uns.get("mean_obs_columns", []))
            sample.adata.uns["mean_obs_columns"] = means + columns
            if hasattr(sample, "clear_multiscale_cache"):
                sample.clear_multiscale_cache()
        self.columns = columns
        if self.component.autocomplete_input is not None:
            self.component.autocomplete_input.completions = self.component.get_colnames()
        if hasattr(self.component, "add_param_form"):
            self.component.add_param_form.refresh_completions()

    def get_widgets(self, doc):
        path = TextInput(title="Single-cell reference (.h5ad path on server)",
                         value=self.reference_path or "")
        load = Button(label="Load reference")
        options = [] if self.reference is None else [str(c) for c in self.reference.obs.columns]
        labels = Select(title="Reference cell-type column", options=options,
                        value=self.label_column if self.label_column in options else (options[0] if options else ""))
        labels.on_change("value", lambda attr, old, new: setattr(self, "label_column", new))
        run = Button(label="Run deconvolution", button_type="success", disabled=not options)
        display = Select(title="Display cell-type fraction", options=self.columns,
                         value=self.columns[0] if self.columns else "", disabled=not self.columns)
        status = LoadingLogger(doc)

        def on_load():
            if status.running:
                return
            reference_path = str(Path(path.value.strip()).expanduser())
            if not path.value.strip() or Path(reference_path).suffix.lower() != ".h5ad":
                status.error("Enter the server path to an annotated .h5ad reference.")
                return
            load.disabled = run.disabled = True

            def done(reference):
                self.reference = reference
                self.reference_path = reference_path
                options = [str(c) for c in reference.obs.columns]
                labels.options = options
                labels.value = options[0] if options else ""
                load.disabled = False
                run.disabled = not options
                status.success(f"Reference loaded: {reference.n_obs:,} cells, {reference.n_vars:,} genes")
                if not options:
                    status.error("Reference has no cell-type annotation columns.")

            def failed(exc):
                load.disabled = False
                run.disabled = self.reference is None or not labels.value
                status.error(exc)

            status.run("Loading reference", lambda: ad.read_h5ad(reference_path), done, failed)

        def on_run():
            if status.running or self.reference is None:
                return
            reference, label_column = self.reference, labels.value
            load.disabled = run.disabled = labels.disabled = True

            def restore():
                load.disabled = run.disabled = labels.disabled = False

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
                       lambda: deconvolve(self.component.sample.adata, reference, label_column), done, failed)

        def show():
            if not display.value:
                return
            if self.component.autocomplete_input is not None:
                self.component.autocomplete_input.value = display.value
            self.component.show_colname(display.value)
            if hasattr(self.component, "loader"):
                self.component.loader.force_reload()

        display.on_change("value", lambda attr, old, new: show())
        load.on_click(on_load)
        run.on_click(on_run)
        return [Div(text="NNLS estimates expression mixture fractions from an annotated single-cell "
                         "reference. Use linear counts in .X and matching gene identifiers. "
                         "Runs on all native observations of this layer; empty observations remain unassigned."),
                path, load, labels, run, status.div, display]
