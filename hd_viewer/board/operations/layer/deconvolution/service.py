"""Validated, cached deconvolution jobs independent of the Bokeh application."""
from contextlib import contextmanager
from dataclasses import dataclass, field
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import threading

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml

from .....plugins.runners import make_runner, RunnerError
from .plugins import ROOT, SDK_DIR, PLUGINS_DIR, discover
from .api.sdk.vhd_deconv import read_table, write_table


@dataclass
class DeconvolutionResult:
    fractions: pd.DataFrame
    genes: list
    abundance: pd.DataFrame = None
    metadata: dict = field(default_factory=dict)

    def __getitem__(self, index):
        return (self.fractions, self.genes)[index]


def load_config(path=ROOT / 'config.yaml'):
    if not Path(path).exists():
        return {}
    return yaml.safe_load(Path(path).read_text()) or {}


def matrix_fingerprint(adata):
    """Canonical sparse representation, including row and gene identities."""
    digest = hashlib.sha256()
    digest.update(json.dumps([list(adata.obs_names), list(adata.var_names)], ensure_ascii=False).encode())
    matrix = sparse.csr_matrix(adata.X, dtype=np.float64, copy=True)
    matrix.sum_duplicates()
    matrix.sort_indices()
    matrix.eliminate_zeros()
    for values in (matrix.data.astype('<f8'), matrix.indices.astype('<i8'), matrix.indptr.astype('<i8')):
        digest.update(values.tobytes())
    return digest.hexdigest()


def _key(spec):
    return hashlib.sha256(json.dumps(spec, sort_keys=True, allow_nan=False).encode()).hexdigest()


def validate_result(result, index, cell_types=None):
    """Reject misaligned or malformed outputs before any viewer/cache mutation."""
    fractions = result.fractions
    if not fractions.index.is_unique or not fractions.index.equals(index):
        raise ValueError('Result observations do not match the target in order.')
    if not fractions.columns.is_unique or not len(fractions.columns):
        raise ValueError('Result cell types must be nonempty and unique.')
    if any(not isinstance(c, str) or not c.strip() or c == 'most_abundant' for c in fractions.columns):
        raise ValueError('Invalid or reserved result cell-type name.')
    if cell_types is not None and set(fractions.columns) != set(cell_types):
        raise ValueError('Result cell types do not match the reference labels.')
    values = fractions.to_numpy(dtype=float)
    missing = np.isnan(values)
    if (missing.any(axis=1) != missing.all(axis=1)).any():
        raise ValueError('Partially missing fraction rows are invalid.')
    if np.isinf(values).any() or (values[~missing] < 0).any():
        raise ValueError('Fractions must be finite and nonnegative (or wholly unassigned).')
    assigned = ~missing.all(axis=1)
    if not np.allclose(values[assigned].sum(axis=1), 1, atol=1e-6, rtol=1e-6):
        raise ValueError('Assigned fractions must sum to one.')
    if result.abundance is not None:
        abundance = result.abundance
        if not abundance.index.equals(index) or not abundance.columns.equals(fractions.columns):
            raise ValueError('Abundance must align with fractions.')
        numbers = abundance.to_numpy(dtype=float)
        if np.isinf(numbers).any() or (numbers[np.isfinite(numbers)] < 0).any():
            raise ValueError('Abundance must be finite and nonnegative.')
        if not np.array_equal(np.isnan(numbers), missing):
            raise ValueError('Abundance assignments must match fractions.')
        if not result.metadata.get('abundance_units'):
            raise ValueError('Abundance units are required.')
        if not np.allclose(numbers[assigned] / numbers[assigned].sum(axis=1, keepdims=True), values[assigned], atol=1e-6):
            raise ValueError('Fractions must be normalized abundance.')
    return result


@contextmanager
def _file_lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


class DeconvolutionService:
    def __init__(self, config=None, *, plugins_dir=PLUGINS_DIR, sdk_dir=SDK_DIR, root=ROOT):
        self.config = load_config() if config is None else dict(config)
        self.plugins = discover(plugins_dir)
        self.runner = make_runner(self.config, sdk_dir=sdk_dir, segmentation_root=root,
                                  env_var='VHD_DECONV_RUNNER', kind='deconvolution')
        self.runtime_available = self.runner.available()
        self.work_dir = Path(self.config.get('work_dir', Path(tempfile.gettempdir()) / f'visium-deconvolution-{os.getuid()}'))
        self._slots = threading.Semaphore(max(1, int(self.config.get('max_parallel', 1))))

    def unavailable_reason(self, name):
        plugin = self.plugins[name]
        if self.runner.name not in plugin.runners:
            return f'{plugin.title} requires Docker or Apptainer; configure a supported deconvolution runner.'
        if not self.runtime_available:
            return f'The configured {self.runner.name} runtime is unavailable; check deconvolution/config.yaml.'
        if name == 'cibersortx' and not all(os.environ.get(k) for k in ('VHD_CIBERSORTX_USERNAME', 'VHD_CIBERSORTX_TOKEN')):
            return 'Set VHD_CIBERSORTX_USERNAME and VHD_CIBERSORTX_TOKEN on the server to enable CIBERSORTx.'
        return ''

    def validate_params(self, name, reference, label_column, params=None):
        plugin = self.plugins[name]
        resolved = plugin.resolve_params(params)
        columns = [label_column] + [resolved[p.name] for p in plugin.params
                                   if p.type == 'obs_column' and resolved[p.name]]
        if not label_column:
            raise ValueError('Choose a reference cell-type column.')
        for column in columns:
            if column not in reference.obs:
                raise ValueError(f'Reference annotation column {column!r} is missing.')
            values = reference.obs[column]
            if values.isna().any() or values.astype(str).str.strip().eq('').any():
                raise ValueError(f'Every reference cell must have a label in {column}.')
        if 'most_abundant' in set(reference.obs[label_column].astype(str)):
            raise ValueError('The cell-type label most_abundant is reserved for the dominant-type readout.')
        for param in plugin.params:
            if param.multiple_subjects and reference.obs[resolved[param.name]].nunique() < 2:
                raise ValueError(f'{param.label} requires at least two reference subjects.')
        if name == 'cibersortx' and resolved['g_min'] > resolved['g_max']:
            raise ValueError('Minimum signature genes must not exceed maximum signature genes.')
        if name == 'music' and resolved['sample_column'] == label_column:
            raise ValueError('Choose a donor/sample column distinct from the cell-type column.')
        return resolved, list(dict.fromkeys(columns))

    def run(self, name, target, reference, label_column, params=None, on_progress=None):
        progress = on_progress or (lambda fraction, message: None)
        reason = self.unavailable_reason(name)
        if reason:
            raise RunnerError(reason)
        plugin = self.plugins[name]
        resolved, columns = self.validate_params(name, reference, label_column, params)
        # Export only expression and required annotations, never stale analysis objects.
        for data in (target, reference):
            if not data.n_obs or not data.n_vars or not data.obs_names.is_unique or not data.var_names.is_unique:
                raise ValueError('Inputs require observations and unique observation/gene identifiers.')
            values = data.X.data if sparse.issparse(data.X) else np.asarray(data.X)
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError('Expression must contain finite, nonnegative linear counts.')
        if len(target.var_names.intersection(reference.var_names)) < 2:
            raise ValueError('Too few shared genes.')
        target = ad.AnnData(target.X, obs=pd.DataFrame(index=target.obs_names.copy()),
                            var=pd.DataFrame(index=target.var_names.copy()))
        reference = ad.AnnData(reference.X, obs=reference.obs[columns].astype(str).copy(),
                               var=pd.DataFrame(index=reference.var_names.copy()))
        reference_spec = dict(matrix=matrix_fingerprint(reference), annotations=reference.obs.to_dict('list'),
                              label_column=label_column, method=name, version=plugin.version,
                              image=plugin.image, runtime=(self.runner.image_uri(plugin) if self.runner.name == 'apptainer'
                                  else self.config.get('local', {}).get('python', {}).get(name, '') if self.runner.name == 'local'
                                  else plugin.image), gpu=self.runner.gpu, params={p: resolved[p] for p in plugin.reference_params})
        reference_key = _key(reference_spec)
        spec = dict(reference=reference_key, target=matrix_fingerprint(target), params=resolved,
                    runner=self.runner.name, gpu=self.runner.gpu)
        key = _key(spec)
        destination = self.work_dir / name / key
        artifacts = self.work_dir / '_references' / name / reference_key
        types = sorted(set(reference.obs[label_column]))
        shared = target.var_names.intersection(reference.var_names, sort=False)
        nonempty = np.asarray(target[:, shared].X.sum(axis=1)).ravel() > 0
        # Reference lock prevents concurrent fits/publication, including other viewer processes.
        progress(0, 'Queued: waiting for a free slot')
        with self._slots, _file_lock(self.work_dir / '_locks' / f'{reference_key}.lock'):
            if (destination / 'result.json').exists():
                result = self._read(destination)
                validate_result(result, target.obs_names, types)
                progress(1, 'Loaded cached result')
                return result
            destination.parent.mkdir(parents=True, exist_ok=True)
            work = Path(tempfile.mkdtemp(prefix=key + '.running-', dir=destination.parent))
            input_dir, output_dir = work / 'input', work / 'output'
            input_dir.mkdir()
            output_dir.mkdir()
            target[nonempty].copy().write_h5ad(input_dir / 'target.h5ad')
            reference.write_h5ad(input_dir / 'reference.h5ad')
            (input_dir / 'params.json').write_text(json.dumps(dict(params=resolved, label_column=label_column)))
            if artifacts.exists():
                shutil.copytree(artifacts, input_dir / 'reference_artifacts')
            try:
                if not nonempty.any():
                    result = DeconvolutionResult(pd.DataFrame(np.nan, index=target.obs_names, columns=types), [])
                    if name == 'cell2location':
                        result.abundance = result.fractions.copy()
                        result.metadata['abundance_units'] = 'expected cells per location (posterior mean)'
                else:
                    if name == 'cibersortx':
                        from .cibersortx import run as run_cibersortx
                        run_cibersortx(plugin, self.runner, input_dir, output_dir, progress, work / 'run.log')
                    else:
                        self.runner.run(plugin, input_dir, output_dir, progress, work / 'run.log')
                    result = self._read(output_dir)
                    validate_result(result, target.obs_names[nonempty], types)
                    result.fractions = result.fractions.reindex(target.obs_names)
                    if result.abundance is not None:
                        result.abundance = result.abundance.reindex(target.obs_names)
                if not isinstance(result.genes, list) or len(set(result.genes)) != len(result.genes) or set(result.genes) - set(target.var_names.intersection(reference.var_names)):
                    raise ValueError('Result genes must be unique shared gene identifiers.')
                result.metadata = {key: result.metadata[key] for key in ('fraction_units', 'abundance_units')
                                   if key in result.metadata}
                result.metadata.update(method=name, version=plugin.version, params=resolved,
                                       reference_id=reference_key, run_id=key, runner=self.runner.name,
                                       log_path=str(destination / 'run.log'))
                validate_result(result, target.obs_names, types)
                self._write(result, output_dir)
                if not (work / 'run.log').exists():
                    (work / 'run.log').write_text('No nonempty target observations.\n')
                shutil.copy(work / 'run.log', output_dir / 'run.log')
                produced = output_dir / 'reference_artifacts'
                if produced.is_dir() and any(produced.iterdir()) and not artifacts.exists():
                    artifacts.parent.mkdir(parents=True, exist_ok=True)
                    staging = Path(tempfile.mkdtemp(dir=artifacts.parent))
                    shutil.copytree(produced, staging, dirs_exist_ok=True)
                    os.replace(staging, artifacts)
                if destination.exists():
                    shutil.rmtree(destination)
                os.replace(output_dir, destination)
                shutil.rmtree(work)
                return result
            except Exception as exc:
                raise RunnerError(f'{plugin.title}: {exc}\nRun files and log: {work}') from exc

    @staticmethod
    def _read(directory):
        metadata = json.loads((directory / 'result.json').read_text())
        abundance = directory / 'abundance.json'
        return DeconvolutionResult(read_table(directory / 'fractions.json'), metadata.pop('genes'),
                                   read_table(abundance) if abundance.exists() else None, metadata)

    @staticmethod
    def _write(result, directory):
        write_table(result.fractions, directory / 'fractions.json')
        if result.abundance is not None:
            write_table(result.abundance, directory / 'abundance.json')
        (directory / 'result.json').write_text(json.dumps(dict(genes=result.genes, **result.metadata), allow_nan=False))
