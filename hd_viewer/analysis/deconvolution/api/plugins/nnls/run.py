from vhd_deconv import start
from nnls import deconvolve

if __name__ == '__main__':
    job = start()
    target, reference = job.read_inputs()
    fractions, genes = deconvolve(target, reference, job.label_column, on_progress=lambda message: job.report(0, message), **job.params)
    job.write_result(fractions.drop(columns='most_abundant'), genes)
