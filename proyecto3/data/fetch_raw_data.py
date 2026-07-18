"""
fetch_raw_data.py — Descarga cruda desde el HEK (fuera del notebook).

Esto NO reemplaza la Etapa 1 del pipeline (esa consulta vive en el notebook
y se ejecuta de la misma forma). Lo hacemos aparte una sola vez porque cada
mes de datos tarda minutos en llegar del servidor del HEK, y no queremos que
"ejecutar el notebook de punta a punta" signifique esperar ~20 minutos cada vez.
El notebook simplemente reutiliza estos CSV crudos si ya existen (igual que
pide la pseudocódigo del enunciado: raw_table.write('raw_data.csv')).
"""

import time
import pandas as pd
from sunpy.net import hek, attrs as a

MESES = [
    ("2014-01-01", "2014-02-01"),
    ("2014-02-01", "2014-03-01"),
    ("2014-03-01", "2014-04-01"),
    ("2014-04-01", "2014-05-01"),
]

client = hek.HEKClient()

dfs_ar = []
dfs_fl = []

for inicio, fin in MESES:
    t0 = time.time()
    print(f"Consultando AR {inicio} .. {fin} ...", flush=True)
    res_ar = client.search(
        a.Time(inicio, fin),
        a.hek.EventType("AR"),
        a.hek.FRM.Name == "NOAA SWPC Observer",
    )
    print(f"  -> {len(res_ar)} filas en {time.time()-t0:.1f} s", flush=True)
    # Convertimos a pandas mes a mes: concatenar Tables de astropy con
    # columnas SkyCoord de distintos meses (distinto obstime) rompe
    # astropy.table.vstack; en pandas es una simple columna de floats.
    dfs_ar.append(res_ar.to_pandas())

    t0 = time.time()
    print(f"Consultando FL {inicio} .. {fin} ...", flush=True)
    res_fl = client.search(
        a.Time(inicio, fin),
        a.hek.EventType("FL"),
    )
    print(f"  -> {len(res_fl)} filas en {time.time()-t0:.1f} s", flush=True)
    dfs_fl.append(res_fl.to_pandas())

df_ar = pd.concat(dfs_ar, ignore_index=True)
df_fl = pd.concat(dfs_fl, ignore_index=True)

df_ar.to_csv("ar_raw.csv", index=False)
df_fl.to_csv("fl_raw.csv", index=False)

print(f"AR: {len(df_ar)} filas -> ar_raw.csv")
print(f"FL: {len(df_fl)} filas -> fl_raw.csv")
print("LISTO")
