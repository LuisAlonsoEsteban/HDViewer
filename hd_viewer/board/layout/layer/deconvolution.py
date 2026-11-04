
from functools import lru_cache
from html import escape
import json
from urllib.parse import quote

import numpy as np
import pandas as pd
from bokeh.layouts import column
from bokeh.models import Button, Checkbox, Div, NumericInput, Select

from .... import settings
from ....settings.config import deconvolution_service_kwargs
from ...operations.layer.deconvolution import DeconvolutionService, DeconvolutionResult, validate_result
from ...utils.loading import LoadingLogger


@lru_cache(maxsize=1)
def _default_service():
    # One queue and runner per viewer process, shared by all layers/sessions.
    return DeconvolutionService(**deconvolution_service_kwargs())


class DeconvolutionManager:
    DECONV_RESULT_PREFIX = settings.DECONV_RESULT_PREFIX

    def __init__(self, component):
        self.component = component
        self.reference = None
        self.reference_path = None
        self.reference_name = None
        self.label_column = ''
        self.columns = []
        self.service = None
        self.work_dir = None
        self.source_bin_size_um = None
        self.method = 'nnls'
        self.parameter_values = {}
        self.current_run = ''

    @property
    def history(self):
        return self.component.sample.adata.uns.get('deconvolution_runs', {})

    @property
    def cache_namespace(self):
        """Native input resolution and layer, independent of display aggregation."""
        sample = self.component.sample
        size = self.source_bin_size_um or getattr(sample, 'bin_size', None)
        resolution = 'native'
        if size is not None:
            size = float(size)
            if np.isfinite(size) and size > 0:
                label = f'{int(size):03d}' if size.is_integer() else f'{size:g}'.replace('.', 'p')
                resolution = f'square_{label}um'
        kind = getattr(self.component, 'sample_type', type(sample).__name__.lower())
        key = getattr(self.component, 'key', kind)
        # Encode identifiers as individual directory names, including dots/slashes.
        parts = [resolution, quote(str(kind), safe='').replace('.', '%2E')]
        if key != kind:
            parts.append(quote(str(key), safe='').replace('.', '%2E'))
        return '/'.join(parts)

    def apply_result(self, result, label_column):
        if not isinstance(result, DeconvolutionResult):
            fractions, genes = result
            result = DeconvolutionResult(fractions.drop(columns='most_abundant', errors='ignore'), genes,
                                         metadata=dict(method='nnls', version='0.1.0', params={}))

        sample = self.component.sample
        validate_result(result, sample.adata.obs_names)
        fractions = result.fractions.copy()
        dominant = pd.Series(np.nan, index=fractions.index, dtype=object)
        assigned = fractions.notna().any(axis=1)
        dominant.loc[assigned] = fractions.loc[assigned].idxmax(axis=1)
        fractions['most_abundant'] = dominant

        # Validate everything before publishing under the sample lock.
        with sample.fetch_lock:
            run = 1
            while str(run) in self.history or any(str(c).startswith(f'{self.DECONV_RESULT_PREFIX}{run}:') for c in sample.get_colnames()):
                run += 1
            prefix = f'{self.DECONV_RESULT_PREFIX}{run}:'
            fraction_columns = [prefix + name for name in result.fractions.columns]
            abundance_columns = []
            for name, original in zip(fraction_columns, result.fractions.columns):
                sample.adata.obs[name] = result.fractions[original].astype(float)
            sample.adata.obs[prefix + 'most_abundant'] = dominant.astype('category')
            if result.abundance is not None:
                for i, original in enumerate(result.abundance.columns):
                    # Separate namespace avoids collisions with biological cell-type names.
                    name = f'{self.DECONV_RESULT_PREFIX}{run}_abundance:{i}'
                    abundance_columns.append(name)
                    sample.adata.obs[name] = result.abundance[original].astype(float)
                sample.adata.obsm[f'deconvolution_{run}_abundance'] = result.abundance.copy()
                sample.adata.obsm['deconvolution_abundance'] = result.abundance.copy()
            else:
                sample.adata.obsm.pop('deconvolution_abundance', None)
            metadata = dict(
                method=result.metadata.get('method', 'nnls'), version=result.metadata.get('version', ''),
                params_json=json.dumps(result.metadata.get('params', {}), sort_keys=True),
                reference_id=result.metadata.get('reference_id', ''), reference_path=self.reference_path or '',
                reference_name=self.reference_name or '', label_column=label_column,
                cell_types=list(result.fractions.columns), genes=result.genes,
                fraction_columns=fraction_columns, abundance_columns=abundance_columns,
                abundance_units=result.metadata.get('abundance_units', ''),
                fraction_units=result.metadata.get('fraction_units', 'relative fraction'),
                log_path=result.metadata.get('log_path', ''), run_id=result.metadata.get('run_id', ''),
                cache_namespace=result.metadata.get('cache_namespace', ''),
            )
            sample.adata.obsm['deconvolution'] = fractions
            sample.adata.obsm[f'deconvolution_{run}'] = fractions.copy()
            sample.adata.uns['deconvolution'] = metadata
            sample.adata.uns.setdefault('deconvolution_runs', {})[str(run)] = metadata.copy()
            means = list(sample.adata.uns.get('mean_obs_columns', []))
            sample.adata.uns['mean_obs_columns'] = list(dict.fromkeys(means + fraction_columns))
            if hasattr(sample, 'clear_multiscale_cache'):
                sample.clear_multiscale_cache()
        self.columns = fraction_columns
        self.current_run = str(run)
        self.component.refresh_completions()

    def get_widgets(self, doc):
        options = [] if self.reference is None else [str(c) for c in self.reference.obs.columns]
        labels = Select(title='Reference cell-type column (.obs)', options=options, disabled=not options,
                        value=self.label_column if self.label_column in options else (options[0] if options else ''))
        labels.on_change('value', lambda attr, old, new: setattr(self, 'label_column', new))
        status = LoadingLogger(doc)
        try:
            self.service = self.service or _default_service()
            plugins = self.service.plugins
            service_error = ''
        except Exception as exc:
            plugins, service_error = {}, str(exc)
        selected = self.method if self.method in plugins else next(iter(plugins), '')
        method = Select(title='Deconvolution method', options=[(p.name, p.title) for p in plugins.values()],
                        value=selected, disabled=not plugins)
        description, availability = Div(), Div()
        parameters = column(sizing_mode='stretch_width')
        parameter_widgets = {}
        run = Button(label='Run deconvolution', button_type='success')
        saved = Select(title='Saved deconvolution run')
        readout = Select(title='Readout', options=[('fractions', 'Relative fractions')], value='fractions')
        display = Select(title='Display cell type', options=[], disabled=True)

        def show():
            if not display.value:
                return
            if self.component.autocomplete_input is not None:
                self.component.autocomplete_input.value = display.value
            self.component.show_colname(display.value)
            if hasattr(self.component, 'loader'):
                self.component.loader.force_reload()

        def update_display():
            record = self.history.get(saved.value, {})
            names = list(record.get('abundance_columns' if readout.value == 'abundance' else 'fraction_columns', []))
            types = list(record.get('cell_types', []))
            display.options = list(zip(names, types))
            display.disabled = not names
            display.value = names[0] if names else ''
            self.columns = names

        def update_readout():
            record = self.history.get(saved.value, {})
            choices = [('fractions', 'Relative fractions')]
            if len(record.get('abundance_columns', [])):
                choices.append(('abundance', 'Abundance: ' + record['abundance_units']))
            readout.options = choices
            if readout.value not in dict(choices):
                readout.value = 'fractions'
            readout.disabled = not record
            update_display()

        def update_runs():
            saved.options = [(key, f"Run {key}: {value['method']}") for key, value in self.history.items()]
            saved.disabled = not self.history
            saved.value = self.current_run if self.current_run in self.history else next(reversed(self.history), '')
            update_readout()

        def raw_params():
            plugin = plugins[method.value]
            values = {}
            for param in plugin.params:
                widget = parameter_widgets[param.name]
                value = widget.active if param.type == 'bool' else widget.value
                if param.type == 'choice':
                    value = next(c for c in param.choices if str(c) == value)
                values[param.name] = value
            return values

        def update_availability():
            if not plugins:
                reason = 'Deconvolution plugins unavailable: ' + service_error
            elif not options:
                reason = ('No single-cell reference was provided for this dataset. Deconvolution is unavailable; '
                          'the other layer tools remain available.' if self.reference is None else
                          'The reference contains counts but no .obs annotation columns. '
                          'Upload a barcode-indexed annotation CSV or Parquet file alongside it in the data loading step.')
            else:
                reason = self.service.unavailable_reason(method.value)
                if not reason:
                    try:
                        self.service.validate_params(method.value, self.reference, labels.value, raw_params())
                    except (ValueError, TypeError) as exc:
                        reason = str(exc)
            availability.text = escape(reason)
            run.disabled = bool(reason) or status.running

        def build_parameters():
            self.method = method.value
            parameter_widgets.clear()
            plugin = plugins.get(method.value)
            if plugin is None:
                update_availability()
                return
            description.text = escape(plugin.description)
            if plugin.url:
                description.text += f' <a href="{escape(plugin.url, quote=True)}" target="_blank">Documentation</a>'
            if plugin.citation:
                description.text += '<br>' + escape(plugin.citation)
            values = self.parameter_values.setdefault(plugin.name, {})
            widgets = []
            for param in plugin.params:
                value = values.get(param.name, param.default)
                if param.type == 'obs_column':
                    widget = Select(title=param.label, options=[('', 'Choose column' if param.required else 'None')] + [(x, x) for x in options],
                                    value=value if value in options else '')
                elif param.type == 'choice':
                    widget = Select(title=param.label, options=[str(x) for x in param.choices], value=str(value))
                elif param.type == 'bool':
                    widget = Checkbox(label=param.label, active=bool(value))
                else:
                    widget = NumericInput(title=param.label, value=value, mode=param.type, low=param.min, high=param.max)
                def change(attr, old, new, name=param.name, store=values):
                    store[name] = new
                    update_availability()
                widget.on_change('active' if param.type == 'bool' else 'value', change)
                parameter_widgets[param.name] = widget
                widgets.append(widget)
            parameters.children = widgets
            update_availability()

        def restore():
            method.disabled = not plugins
            labels.disabled = not options
            for widget in parameter_widgets.values():
                widget.disabled = False
            update_availability()

        def on_run():
            if status.running or run.disabled:
                return
            name, label_column = method.value, labels.value
            try:
                params, _ = self.service.validate_params(name, self.reference, label_column, raw_params())
            except Exception as exc:
                status.error(exc)
                return
            # Snapshot counts/annotations while on the document thread. Later UI edits
            # cannot change which input data a background run consumes.
            sample = self.component.sample
            with sample.fetch_lock:
                import anndata as ad
                target = ad.AnnData(sample.adata.X.copy(), obs=pd.DataFrame(index=sample.adata.obs_names.copy()),
                                    var=pd.DataFrame(index=sample.adata.var_names.copy()))
                if name == 'spacexr' and 'spatial' in sample.adata.obsm:
                    target.obsm['spatial'] = sample.adata.obsm['spatial'].copy()
            reference = self.reference
            work_dir = self.work_dir
            cache_namespace = self.cache_namespace
            run.disabled = method.disabled = labels.disabled = True
            for widget in parameter_widgets.values():
                widget.disabled = True

            def done(result):
                try:
                    self.apply_result(result, label_column)
                    update_runs()
                    status.success(f'Deconvolved {len(result.fractions):,} observations into {len(result.fractions.columns)} cell types')
                except Exception as exc:
                    status.error(exc)
                finally:
                    restore()

            def failed(exc):
                status.error(exc)
                restore()

            status.run('Deconvolution (' + plugins[name].title + ')',
                       lambda: self.service.run(name, target, reference, label_column, params,
                           work_dir=work_dir, cache_namespace=cache_namespace,
                           on_progress=lambda fraction, message: status.log(message)), done, failed)

        labels.on_change('value', lambda attr, old, new: update_availability())
        method.on_change('value', lambda attr, old, new: build_parameters())
        saved.on_change('value', lambda attr, old, new: update_readout())
        readout.on_change('value', lambda attr, old, new: update_display())
        display.on_change('value', lambda attr, old, new: show())
        run.on_click(on_run)
        build_parameters()
        update_runs()
        if self.reference is not None and options:
            status.set_message(escape(f'Reference ready: {self.reference_name or "loaded reference"} '
                                     f'({self.reference.n_obs:,} cells, {self.reference.n_vars:,} genes). '
                                     f'Runner: {self.service.runner.name if self.service else "unavailable"}.'))
        else:
            status.set_message(availability.text)
        return [Div(text='Run a method on all native observations of this layer. Use linear counts and matching gene identifiers. '
                         'Empty observations remain unassigned.'), labels, method, description, parameters, availability,
                run, status.div, saved, readout, display]
