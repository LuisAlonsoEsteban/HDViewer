
from html import escape

import numpy as np
from bokeh.models import Button, Checkbox, CheckboxGroup, Div, MultiChoice, Select

from ...operations.layer.metadata_comparison import compare_metadata, compare_metadata_to_rest, rest_mask, numeric_columns, build_metadata_plot
from ...utils.loading import LoadingLogger


class MetadataComparisonManager:
    def __init__(self, component):
        self.component = component
        self._selected = {'groups': set(), 'clusters': set()}
        self._checkboxes = {}
        self._generation = 0
        self._context = None
        self.result = None
        self._layout_built = False
        self._refreshing = False
        component.groups_manager.add_change_listener(self._groups_changed)

    @property
    def sample(self):
        return self.component.sample

    @property
    def selectors(self):
        return self.component.de_manager

    def _presets(self):
        records = dict(self.sample.adata.uns.get('deconvolution_runs', {}))
        if not records and 'deconvolution' in self.sample.adata.uns:
            records['latest'] = self.sample.adata.uns['deconvolution']
        specs = {}
        numeric = set(numeric_columns(self.sample.adata.obs))
        ordered = sorted(records.items(), key=lambda item: (1, int(item[0])) if str(item[0]).isdigit() else (0, str(item[0])))
        for run, record in ordered:
            for kind, title, unit in (('fraction', 'fractions', 'proportion'),
                                      ('abundance', 'abundance', str(record.get('abundance_units', 'abundance')))):
                pairs = [(str(col), str(label)) for col, label in zip(record.get(kind + '_columns', []), record.get('cell_types', []))
                         if col in numeric]
                if pairs:
                    specs[f'{run}:{kind}'] = dict(title=f"Run {run}: {record.get('method', 'deconvolution')} {title}",
                                                  labels=dict(pairs), unit=unit, fractions=kind == 'fraction')
        return specs

    def refresh_columns(self):
        if not self._layout_built:
            return
        specs = self._presets()
        previous = self.columns.value.copy()
        current = self.preset.value
        options = [('numeric', 'Numeric metadata')] + [(key, spec['title']) for key, spec in specs.items()]
        if current not in dict(options):
            current = 'numeric'
        if not previous and specs and current == 'numeric':
            current = next(reversed(specs))
            # Prefer fractions when both readouts exist for the latest run.
            fraction = current.rsplit(':', 1)[0] + ':fraction'
            current = fraction if fraction in specs else current
        self._refreshing = True
        self.preset.options = options
        self.preset.value = current
        self._refreshing = False
        self._set_columns(previous)
        if not self._context_matches():
            self._invalidate('Metadata comparison cleared because its columns or observations changed.')

    def _set_columns(self, previous=None):
        specs = self._presets()
        spec = specs.get(self.preset.value)
        pairs = list(spec['labels'].items()) if spec else [(name, name) for name in numeric_columns(self.sample.adata.obs)]
        allowed = dict(pairs)
        selected = [name for name in (self.columns.value if previous is None else previous) if name in allowed]
        self.columns.options = pairs
        self.columns.value = selected if selected else list(allowed) if spec else []
        self.columns.disabled = not pairs
        self.hint.text = ('Uses all native observations, independently of the viewport or displayed bin aggregation. '
                          'Overlapping observations contribute to each group. Missing and infinite values are omitted per column; '
                          'bars show descriptive summaries, without significance tests. '
                          'Select The rest for one plot per selected group/cluster against the observations outside all defined custom groups '
                          'and outside that group/cluster.')

    def _refresh_sets(self):
        if not self._layout_built:
            return
        self._refreshing = True
        for source, checkbox in self._checkboxes.items():
            names = self.selectors._de_options(source)
            self._selected[source].intersection_update(names)
            checkbox.labels = names
            checkbox.active = [i for i, name in enumerate(names) if name in self._selected[source]]
        self._refreshing = False

    def _context_matches(self):
        context = self._context
        if context is None:
            return True
        adata = self.sample.adata
        if adata is not context['adata'] or not adata.obs_names.equals(context['names']):
            return False
        for source, name in context['chosen']:
            if name not in self.selectors._de_options(source):
                return False
            mask = self.selectors._de_masks(source, [name]).get(name)
            if mask is None or not np.array_equal(mask, context['masks'][name]):
                return False
        if context.get('defined_masks') is not None:
            names = self.selectors._de_options('groups')
            if set(names) != set(context['defined_masks']):
                return False
            current = self.selectors._de_masks('groups', names)
            if set(current) != set(context['defined_masks']) or any(
                    not np.array_equal(current[name], mask) for name, mask in context['defined_masks'].items()):
                return False
        return all(col in adata.obs for col in context['columns'])

    def _invalidate(self, message):
        self._generation += 1
        self._context = None
        self.result = None
        layout = getattr(self.component, '_figure_de_layout', None)
        if layout is not None and getattr(layout, '_comparison_owner', None) is self:
            layout.children = []
            layout._comparison_owner = None
        if self._layout_built and not self.logger.running:
            self.logger.set_message(message)

    def _groups_changed(self):
        self._refresh_sets()
        if not self._context_matches():
            self._invalidate('Metadata comparison cleared because the selected groups changed.')

    def refresh_clusters(self):
        self._selected['clusters'].clear()
        if self._context and any(source == 'clusters' for source, name in self._context['chosen']):
            self._invalidate('Metadata comparison cleared because clustering changed.')
        self._refresh_sets()

    def get_widgets(self, doc):
        self.logger = LoadingLogger(doc)
        self.preset = Select(title='Metadata set', value='numeric', options=[('numeric', 'Numeric metadata')])
        self.columns = MultiChoice(title='Metadata columns', value=[], options=[])
        self.statistic = Select(title='Summary', options=[('mean', 'Mean'), ('median', 'Median')], value='mean')
        self.include_rest = Checkbox(label='The rest (outside all defined custom groups and the compared set)', active=False)
        self.hint = Div()
        self._checkboxes = {source: CheckboxGroup(labels=[], active=[]) for source in ('groups', 'clusters')}
        self.run = Button(label='Compare metadata', button_type='primary')
        self._layout_built = True
        self._refresh_sets()
        self.refresh_columns()
        for source, checkbox in self._checkboxes.items():
            def changed(attr, old, new, source=source, checkbox=checkbox):
                if not self._refreshing:
                    self._selected[source] = {checkbox.labels[i] for i in new}
            checkbox.on_change('active', changed)
        self.preset.on_change('value', lambda attr, old, new: None if self._refreshing else self._set_columns([]))
        self.run.on_click(lambda: self._run(doc))
        return [Div(text='Compare selected groups/clusters, or compare each selected set against The rest.'),
                Div(text='<b>Groups</b>'), self._checkboxes['groups'], Div(text='<b>Clusters</b>'),
                self._checkboxes['clusters'], self.include_rest, self.preset, self.columns, self.statistic,
                self.hint, self.run, self.logger.div]

    def _run(self, doc):
        if self.logger.running:
            return
        chosen = [(source, name) for source in ('groups', 'clusters')
                  for name in self.selectors._de_options(source) if name in self._selected[source]]
        include_rest = self.include_rest.active
        if len(chosen) < (1 if include_rest else 2) or not self.columns.value:
            requirement = 'at least one group or cluster with The rest' if include_rest else 'at least two groups or clusters'
            self.logger.error(f'Select {requirement} and one or more numeric metadata columns.')
            return
        masks, colors = {}, {}
        for source, name in chosen:
            masks.update(self.selectors._de_masks(source, [name]))
            colors.update(self.selectors._de_colors(source))
        if len(masks) != len(chosen):
            self.logger.error('Some selected groups no longer have membership data. Refresh the group selection.')
            return
        columns = self.columns.value.copy()
        statistic = self.statistic.value
        spec = self._presets().get(self.preset.value, {})
        with self.sample.fetch_lock:
            adata = self.sample.adata
            masks = {name: mask.copy() for name, mask in masks.items()}
            defined_masks = None
            union = np.logical_or.reduce(list(masks.values()))
            if include_rest:
                defined_names = self.selectors._de_options('groups')
                defined_masks = {name: mask.copy() for name, mask in self.selectors._de_masks('groups', defined_names).items()}
                if set(defined_names) != set(defined_masks):
                    self.logger.error('Some defined groups have no membership data. Refresh the group definitions.')
                    return
                union |= rest_mask(defined_masks, adata.n_obs)
            frame = adata.obs.loc[union, columns].copy()
            names = adata.obs_names.copy()
        self._context = dict(adata=adata, names=names, masks=masks, chosen=chosen, columns=columns, defined_masks=defined_masks)
        generation = self._generation
        self.run.disabled = True

        def done(result):
            try:
                if generation != self._generation or not self._context_matches():
                    self.logger.set_message('Metadata results discarded because the comparison data changed.')
                    return
                if not adata.obs.loc[frame.index, columns].equals(frame):
                    self._invalidate('Metadata results discarded because column values changed.')
                    return
                layout = getattr(self.component, '_figure_de_layout', None)
                if layout is None:
                    raise ValueError('The figure results panel is unavailable.')
                comparisons = result if include_rest else {'Selected groups': result}
                children = [Div(text=f'<b>Metadata comparison · {escape(self.component.name)}</b>')]
                for name, comparison in comparisons.items():
                    pair_colors = colors.copy()
                    if include_rest:
                        pair_colors[comparison.groups[-1]] = '#888888'
                    chart = build_metadata_plot(comparison, spec.get('labels'), pair_colors,
                                                spec.get('unit', 'value'), spec.get('fractions', False))
                    if include_rest:
                        chart.title.text = f'{name} vs {comparison.groups[-1]}'
                    children.append(chart)
                layout._comparison_owner = self
                layout.children = children
                self.result = result
                self.logger.success('Metadata comparison ready below the map.')
            except Exception as exc:
                self.logger.error(exc)
            finally:
                self.run.disabled = False

        def failed(exc):
            self.run.disabled = False
            self.logger.error(exc)

        def work():
            selected = {name: mask[union] for name, mask in masks.items()}
            if include_rest:
                defined = {name: mask[union] for name, mask in defined_masks.items()}
                return compare_metadata_to_rest(frame, selected, defined, columns, statistic)
            return compare_metadata(frame, selected, columns, statistic)

        self.logger.run('Metadata comparison', work, done, failed)

    def dispose(self):
        self.component.groups_manager.remove_change_listener(self._groups_changed)
        self._invalidate('Layer removed.')
