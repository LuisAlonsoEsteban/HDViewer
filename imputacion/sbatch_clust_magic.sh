#!/bin/bash
#SBATCH --job-name=t_clust_magic
#SBATCH --output=logs/t_clustering_dom_magic_%j.out
#SBATCH --error=logs/t_clustering_dom_magic_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=15G
#SBATCH --time=02:00:00

# Evalua UN SOLO metodo (--metodo magic). Cada job escribe en su propio CSV
# benchmark_*_dominios_magic.csv: con un CSV compartido los cuatro jobs en
# paralelo se pisaban el fichero y solo sobrevivia el ultimo en terminar.

METODO=magic
PYTHON_BIN=/scratch/tcondegerma/conda-envs/tfm-impute/bin/python
cd /scratch/tcondegerma/VisiumVisualizer/imputacion

TIMEFILE=logs/recursos_${METODO}_${SLURM_JOB_ID}.time

echo "Python en uso: $PYTHON_BIN"
$PYTHON_BIN -c "import anndata, numpy, scipy, scanpy, sklearn, psutil; print('imports OK')" || { echo "FALLO IMPORTS"; exit 1; }

echo "==================================================================="
echo "JOB: $SLURM_JOB_ID   NODO: $(hostname)   FECHA: $(date)"
echo "SCRIPT: evaluar_clustering_metodos.py --metodo $METODO --datasets nb_s00_dom_libre nb_s10_dom_libre"
echo "==================================================================="

# -o a fichero propio: en el .err el resumen quedaba sepultado bajo miles de
# lineas de tqdm y avisos de TensorFlow. El coste POR METODO (t_wall/RSS de la
# llamada de imputacion sola) lo mide el propio script Python y va al CSV;
# esto de aqui es el coste del job entero (incluye imports, carga y Leiden).
/usr/bin/time -v -o "$TIMEFILE" $PYTHON_BIN evaluar_clustering_metodos.py --metodo $METODO --datasets nb_s00_dom_libre nb_s10_dom_libre
RC=$?

echo "==================================================================="
echo "Fin: $(date)   codigo de salida: $RC"
echo "RECURSOS DEL JOB COMPLETO (/usr/bin/time -v):"
grep -E "Elapsed \(wall|User time|System time|Maximum resident|Percent of CPU" "$TIMEFILE"
echo "==================================================================="
# sacct NO se consulta aqui: mientras el job corre el step sigue en RUNNING y
# MaxRSS/TotalCPU salen vacios (era el motivo de las columnas en blanco).
# Cuando terminen todos:  ./recursos_jobs.sh <jobid> [<jobid> ...]
echo "Contabilidad de SLURM:  ./recursos_jobs.sh $SLURM_JOB_ID"
exit $RC
