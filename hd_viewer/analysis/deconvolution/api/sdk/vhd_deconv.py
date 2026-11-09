"""File contract for isolated deconvolution plugins (no viewer dependency)."""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def read_table(path):
    # JSON preserves numeric-looking barcodes and cell-type labels exactly.
    data = json.loads(Path(path).read_text())
    return pd.DataFrame(data['data'], index=pd.Index(data['index']), columns=data['columns'], dtype=float)


def write_table(frame, path):
    data = frame.to_numpy(dtype=float)
    if np.isinf(data).any():
        raise ValueError('Infinite result values')
    values = [[None if np.isnan(x) else float(x) for x in row] for row in data]
    Path(path).write_text(json.dumps(dict(index=frame.index.tolist(), columns=frame.columns.tolist(), data=values), allow_nan=False))


class Job:
    def __init__(self, input_dir, output_dir):
        self.input_dir, self.output_dir = Path(input_dir), Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        spec = json.loads((self.input_dir / 'params.json').read_text())
        self.params = spec['params']
        self.label_column = spec['label_column']
        self.reference_cache = self.input_dir / 'reference_artifacts'
        self.reference_output = self.output_dir / 'reference_artifacts'
        self.reference_output.mkdir(exist_ok=True)

    def read_inputs(self):
        return ad.read_h5ad(self.input_dir / 'target.h5ad'), ad.read_h5ad(self.input_dir / 'reference.h5ad')

    def report(self, fraction, message=''):
        print(f'##progress {float(fraction):.4f} {message}', flush=True)

    def write_result(self, fractions, genes, *, abundance=None, abundance_units=''):
        write_table(fractions, self.output_dir / 'fractions.json')
        if abundance is not None:
            write_table(abundance, self.output_dir / 'abundance.json')
        (self.output_dir / 'result.json').write_text(json.dumps(dict(
            genes=list(genes), fraction_units='relative fraction', abundance_units=abundance_units,
        )))
        self.report(1, 'Result written')


def start(argv=None):
    parser = argparse.ArgumentParser(description='Visium deconvolution plugin')
    parser.add_argument('--input', default='/data/input')
    parser.add_argument('--output', default='/data/output')
    args = parser.parse_args(argv)
    
    return Job(args.input, args.output)
