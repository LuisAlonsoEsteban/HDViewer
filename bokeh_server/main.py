
from bokeh.plotting import curdoc
from bokeh.layouts import column

import pandas as pd  # <-- Importamos pandas aquí

# =====================================================================
# EXPERIMENTO DE CONTROL: Imprimir límites del parquet en la terminal
# =====================================================================
try:
    path_parquet = "/scratch/tcondegerma/data/Visium_HD_Mouse_Small_Intestine/binned_outputs/square_016um/spatial/tissue_positions.parquet"
    df_pos = pd.read_parquet(path_parquet)

    print("\n" + "="*50)
    print("[DEBUG EXPERIMENTO] RANGO REAL EN PARQUET:")
    print(f"   -> X (pxl_col): min={df_pos['pxl_col_in_fullres'].min()}, max={df_pos['pxl_col_in_fullres'].max()}")
    print(f"   -> Y (pxl_row): min={df_pos['pxl_row_in_fullres'].min()}, max={df_pos['pxl_row_in_fullres'].max()}")
    print("="*50 + "\n")
except Exception as e:
    print(f"\n No se pudo leer el parquet del experimento: {e}\n")
# =====================================================================


from .panels.startup import create_startup


doc = curdoc()
create_startup(doc)

placeholders = {
    "metadata": column(),
    "embedding_panel": column(),
    "spots_widgets": column(),
    "bins_widgets": column(),
    "cells_widgets": column(),
    "plot_column": column(sizing_mode='stretch_width'),
    "thumbnail_plot": column(sizing_mode='stretch_width'),
}

# Add them to the document with names
for name, comp in placeholders.items():
    comp.name = name
    doc.add_root(comp)
