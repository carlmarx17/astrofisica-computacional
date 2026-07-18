"""
generar_notebook.py — Construye el cuaderno Jupyter de forma programática.

Igual que en el proyecto 1, escribimos el cuaderno como un script de Python
en lugar de editar el JSON del .ipynb a mano: el código fuente es legible,
el diff de git es útil, y regenerar el cuaderno es un solo comando:

    python generar_notebook.py
"""

import nbformat as nbf

cuaderno = nbf.v4.new_notebook()
celdas = []

def md(texto):
    celdas.append(nbf.v4.new_markdown_cell(texto))

def code(texto):
    celdas.append(nbf.v4.new_code_cell(texto))


# ═════════════════════════════════════════════════════════════════════════
# Portada
# ═════════════════════════════════════════════════════════════════════════

md(r"""# 🔆 Regiones Activas Solares — De la Consulta al Dashboard
*Astrofísica Computacional · Proyecto 3 · Opción A (HEK / SunPy)*

---

Una región activa (AR, *Active Region*) es una zona de la fotosfera solar donde
un tubo de flujo magnético emerge desde la zona de convección y forma manchas
solares. Su topología magnética — clasificada por NOAA en el esquema de Hale
como α, β, γ, δ y combinaciones (β-γ, β-δ, β-γ-δ) — determina qué tan
probable es que la región produzca una fulguración.

En este proyecto construimos un pipeline completo:

1. **Consultamos** el *Heliophysics Event Knowledgebase* (HEK) vía SunPy/Fido
   para obtener regiones activas (NOAA SWPC) y fulguraciones (GOES) reales.
2. **Diseñamos y poblamos** una base de datos SQLite normalizada, con una fila
   por región activa y sus parámetros físicos derivados.
3. **Consultamos la base con SQL** (agregaciones, filtros, `GROUP BY`) desde
   pandas.
4. **Construimos un dashboard estadístico** de 6 paneles con interpretación
   física de cada uno.

**Ventana de datos:** 2014-01-01 a 2014-05-01 (parte del máximo del ciclo
solar 24 — buena estadística de fulguraciones, incluida la fulguración
X4.9 del 25 de febrero de 2014).
""")


# ═════════════════════════════════════════════════════════════════════════
# Importaciones
# ═════════════════════════════════════════════════════════════════════════

code(r"""%matplotlib inline
import os
import sqlite3
import datetime
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

plt.rcParams.update({'font.size': 11, 'figure.dpi': 110})
sns.set_style('whitegrid')
warnings.filterwarnings('ignore', category=FutureWarning)

RUTA_DATOS = 'data'
RUTA_DB    = 'astro_project.db'
os.makedirs(RUTA_DATOS, exist_ok=True)

FECHA_INICIO = '2014-01-01'
FECHA_FIN    = '2014-05-01'   # 4 meses, cerca del máximo del ciclo solar 24

# Conversión de área: NOAA/SWPC reporta el área de manchas en "millonésimas
# de hemisferio solar" (MSH). 1 hemisferio solar = 2π R_sol².
RADIO_SOL_KM = 696_000.0
MM2_POR_MSH  = (2 * np.pi * RADIO_SOL_KM**2) / 1e6 * 1e-6   # Mm² por 1 MSH
print(f"1 MSH equivale a {MM2_POR_MSH:.3f} Mm^2")
""")


# ═════════════════════════════════════════════════════════════════════════
# ETAPA 1 — Adquisición de datos (HEK)
# ═════════════════════════════════════════════════════════════════════════

md(r"""## Etapa 1 · Consultar el HEK

El HEK cataloga eventos solares detectados por distintos módulos ("features
recognition modules", FRM). Para el mismo tipo de evento (`AR` = región
activa) coexisten varios catálogos automáticos:

| `frm_name` | Qué es | ¿Sirve para este proyecto? |
|---|---|---|
| `SPoCA` | Detección automática por color en EUV | No da número NOAA ni clase magnética |
| `HMI SHARP` | Recortes automáticos de magnetogramas | No da clase magnética de Hale |
| **`NOAA SWPC Observer`** | Reporte diario oficial de NOAA/SWPC | **Sí**: número NOAA, clase de Hale (α/β/γ/δ), número de manchas, área |

Por eso filtramos explícitamente `frm_name == 'NOAA SWPC Observer'`.

Para las fulguraciones asociadas (conteo y flujo pico de rayos X) consultamos
por separado el tipo de evento `FL`, que trae `ar_noaanum` (para unir con las
regiones activas) y `fl_peakflux` en W/m² (el flujo GOES estándar).

**Nota de reproducibilidad:** cada mes de consulta al HEK tarda del orden de
minutos (el servidor es lento, no nuestro código). Para que el cuaderno se
pueda ejecutar de punta a punta en segundos y no en ~15 minutos, cacheamos la
respuesta cruda en `data/ar_raw.csv` y `data/fl_raw.csv` — exactamente lo que
sugiere el pseudocódigo del enunciado (`raw_table.write('raw_data.csv')`). Si
borras esos archivos, la celda vuelve a consultar el HEK en vivo con el mismo
código.
""")

code(r"""class QueryError(Exception):
    '''Se lanza cuando el HEK no devuelve datos utilizables.'''
    pass


def consultar_hek_por_meses(tipo_evento, filtro_frm=None, fecha_ini=FECHA_INICIO,
                             fecha_fin=FECHA_FIN):
    '''Consulta el HEK mes a mes (el servidor es más confiable con
    ventanas cortas) y concatena los resultados en un solo DataFrame.'''
    from sunpy.net import hek, attrs as a
    from astropy.table import vstack

    cliente = hek.HEKClient()
    meses = pd.date_range(fecha_ini, fecha_fin, freq='MS')
    tablas = []

    for inicio, fin in zip(meses[:-1], meses[1:]):
        args = [a.Time(inicio.date().isoformat(), fin.date().isoformat()),
                a.hek.EventType(tipo_evento)]
        if filtro_frm is not None:
            args.append(a.hek.FRM.Name == filtro_frm)
        resultado = cliente.search(*args)
        print(f"  {tipo_evento} {inicio.date()} .. {fin.date()}: {len(resultado)} filas")
        if len(resultado) > 0:
            tablas.append(resultado)

    if not tablas:
        raise QueryError(f"El HEK no devolvió eventos de tipo {tipo_evento!r}")

    return vstack(tablas).to_pandas()


ruta_ar_cache = os.path.join(RUTA_DATOS, 'ar_raw.csv')
ruta_fl_cache = os.path.join(RUTA_DATOS, 'fl_raw.csv')

if os.path.exists(ruta_ar_cache) and os.path.exists(ruta_fl_cache):
    print("Usando caché local (data/ar_raw.csv, data/fl_raw.csv).")
    df_ar_crudo = pd.read_csv(ruta_ar_cache)
    df_fl_crudo = pd.read_csv(ruta_fl_cache)
else:
    print("No hay caché — consultando el HEK en vivo (varios minutos)...")
    df_ar_crudo = consultar_hek_por_meses('AR', filtro_frm='NOAA SWPC Observer')
    df_fl_crudo = consultar_hek_por_meses('FL')
    df_ar_crudo.to_csv(ruta_ar_cache, index=False)
    df_fl_crudo.to_csv(ruta_fl_cache, index=False)

if len(df_ar_crudo) == 0:
    raise QueryError('No se recuperaron regiones activas — revisa los parámetros')

columnas_esperadas_ar = {'ar_noaanum', 'ar_noaaclass', 'hgs_x', 'hgs_y',
                          'area_atdiskcenter', 'ar_numspots', 'event_starttime'}
faltantes = columnas_esperadas_ar - set(df_ar_crudo.columns)
assert not faltantes, f"Faltan columnas esperadas en el HEK: {faltantes}"

n_antes = len(df_ar_crudo)

# Cortes de calidad: sin número NOAA, sin clase magnética o sin coordenadas
# heliográficas, el registro no sirve para el análisis físico.
df_ar_limpio = df_ar_crudo.dropna(subset=['ar_noaanum', 'ar_noaaclass', 'hgs_x', 'hgs_y'])
df_ar_limpio = df_ar_limpio[df_ar_limpio['ar_noaanum'] > 0]

n_despues = len(df_ar_limpio)
print(f"Regiones activas (filas diarias): {n_antes} -> {n_despues} tras cortes de calidad")
print(f"Número de regiones activas NOAA distintas: {df_ar_limpio['ar_noaanum'].nunique()}")
""")

md(r"""**Preguntas del enunciado:**

- *¿Cuántos objetos devolvió la consulta antes y después de los cortes de
  calidad?* — Ver la salida de la celda anterior. Las filas son observaciones
  **diarias** de cada AR (NOAA la reporta una vez por día mientras cruza el
  disco visible), no una fila por región; por eso agregamos por
  `ar_noaanum` en la siguiente sección.
- *¿Sesgos de selección conocidos?* Sí, dos importantes: (1) NOAA/SWPC solo
  reporta regiones en el **hemisferio visible desde la Tierra** — una región
  puede emerger, evolucionar y decaer en la cara oculta sin ser jamás
  catalogada; (2) el catálogo favorece regiones con manchas suficientemente
  grandes para ser reconocidas de forma consistente; regiones muy pequeñas o
  efímeras (clase magnética simple, de vida corta) están subrepresentadas.
- *¿Qué ventana cubre la muestra?* 2014-01-01 a 2014-05-01, cuatro meses
  cerca del máximo del ciclo solar 24, limitados al disco visible (no hay
  restricción de "cielo" como en un survey estelar — el "campo de visión" es
  literalmente la cara visible del Sol, que rota con período sinódico de
  ~27.3 días).
""")


# ═════════════════════════════════════════════════════════════════════════
# Agregación: una fila por región activa
# ═════════════════════════════════════════════════════════════════════════

md(r"""## De observaciones diarias a un catálogo por región

Cada NOAA AR aparece muchas veces en `df_ar_limpio` (una fila por día que
NOAA la observó). Para el análisis estadístico necesitamos **una fila por
región activa**, con parámetros representativos de toda su vida visible:

- **Área [Mm²]:** el máximo observado (el área pico durante su tránsito,
  convertida de MSH a Mm²).
- **Coordenadas heliográficas de Stonyhurst [°]:** las del instante más
  cercano al meridiano central (`|hgs_x|` mínimo) — la medición geométrica
  menos afectada por el escorzo de proyección cerca del limbo.
- **Clase magnética de Hale:** la moda (clasificación más frecuente)
  observada durante su tránsito — la clasificación puede evolucionar
  día a día conforme la región se desarrolla.
- **Número de manchas:** el máximo observado.
- **Número de fulguraciones y flujo GOES pico [W/m²]:** se obtienen uniendo
  con la tabla de fulguraciones (`FL`) por `ar_noaanum`. Una región sin
  fulguraciones asociadas queda con conteo 0 y flujo pico nulo (`NaN`) — esto
  es información física real, no un dato faltante por error.
""")

code(r"""df_ar_limpio = df_ar_limpio.copy()
df_ar_limpio['area_mm2'] = df_ar_limpio['area_atdiskcenter'] * MM2_POR_MSH
df_ar_limpio['ar_noaanum'] = df_ar_limpio['ar_noaanum'].astype(int)

# ── Coordenadas representativas: la observación más cercana al meridiano central ──
idx_meridiano = df_ar_limpio.groupby('ar_noaanum')['hgs_x'].apply(lambda s: s.abs().idxmin())
coords_centrales = df_ar_limpio.loc[idx_meridiano, ['ar_noaanum', 'hgs_x', 'hgs_y', 'event_starttime']]
coords_centrales = coords_centrales.rename(columns={'hgs_x': 'hgs_lon_deg',
                                                     'hgs_y': 'hgs_lat_deg',
                                                     'event_starttime': 'fecha_meridiano'})

# ── Área pico y número de manchas pico ──────────────────────────────────────
picos = df_ar_limpio.groupby('ar_noaanum').agg(
    area_mm2=('area_mm2', 'max'),
    num_manchas=('ar_numspots', 'max'),
    primera_observacion=('event_starttime', 'min'),
).reset_index()

# ── Clase magnética dominante (moda) ────────────────────────────────────────
clase_dominante = (df_ar_limpio.groupby('ar_noaanum')['ar_noaaclass']
                   .agg(lambda s: s.value_counts().idxmax())
                   .reset_index(name='clase_magnetica'))

df_regiones = picos.merge(coords_centrales, on='ar_noaanum').merge(clase_dominante, on='ar_noaanum')

# ── Fulguraciones asociadas: unir con la tabla FL por número NOAA ───────────
df_fl_limpio = df_fl_crudo.dropna(subset=['ar_noaanum', 'fl_peakflux'])
df_fl_limpio = df_fl_limpio[df_fl_limpio['ar_noaanum'] > 0].copy()
df_fl_limpio['ar_noaanum'] = df_fl_limpio['ar_noaanum'].astype(int)

resumen_flares = df_fl_limpio.groupby('ar_noaanum').agg(
    num_fulguraciones=('fl_peakflux', 'count'),
    flujo_pico_goes=('fl_peakflux', 'max'),
).reset_index()

df_regiones = df_regiones.merge(resumen_flares, on='ar_noaanum', how='left')
df_regiones['num_fulguraciones'] = df_regiones['num_fulguraciones'].fillna(0).astype(int)

# ── Rotación de Carrington de la primera observación (para la serie temporal) ──
from sunpy.coordinates.sun import carrington_rotation_number
df_regiones['rotacion_carrington'] = df_regiones['primera_observacion'].apply(
    lambda t: int(carrington_rotation_number(t))
)

print(f"Catálogo final: {len(df_regiones)} regiones activas distintas")
print(f"Con al menos una fulguración asociada: {(df_regiones['num_fulguraciones'] > 0).sum()}")
df_regiones.head()
""")


# ═════════════════════════════════════════════════════════════════════════
# ETAPA 2 — Diseño y población de la base SQLite
# ═════════════════════════════════════════════════════════════════════════

md(r"""## Etapa 2 · Base de datos SQLite

Diseñamos una tabla `active_regions` con un renglón por región activa y una
tabla `metadata` que documenta de dónde salió el catálogo (fuente, ventana de
tiempo, cortes de calidad aplicados) — así la base es auto-descriptiva y
trazable, tal como sugiere el enunciado.
""")

code(r"""conn = sqlite3.connect(RUTA_DB)
cur = conn.cursor()

cur.execute('''
CREATE TABLE IF NOT EXISTS active_regions (
    id                   INTEGER PRIMARY KEY AUTOINCREMENT,
    noaa_ar              INTEGER NOT NULL UNIQUE,
    hgs_lon_deg          REAL,
    hgs_lat_deg          REAL,
    area_mm2             REAL,
    clase_magnetica      TEXT,
    num_manchas          INTEGER,
    num_fulguraciones    INTEGER,
    flujo_pico_goes      REAL,
    rotacion_carrington  INTEGER,
    fecha_observacion    TEXT,
    retrieved_at         TEXT
)
''')

cur.execute('''
CREATE TABLE IF NOT EXISTS metadata (
    key   TEXT PRIMARY KEY,
    value TEXT
)
''')

marca_tiempo = datetime.datetime.now(datetime.timezone.utc).isoformat()

registros = [
    (
        int(fila.ar_noaanum), float(fila.hgs_lon_deg), float(fila.hgs_lat_deg),
        float(fila.area_mm2), str(fila.clase_magnetica), int(fila.num_manchas),
        int(fila.num_fulguraciones),
        None if pd.isna(fila.flujo_pico_goes) else float(fila.flujo_pico_goes),
        int(fila.rotacion_carrington), str(fila.primera_observacion), marca_tiempo,
    )
    for fila in df_regiones.itertuples()
]

cur.executemany('''
    INSERT OR REPLACE INTO active_regions
        (noaa_ar, hgs_lon_deg, hgs_lat_deg, area_mm2, clase_magnetica,
         num_manchas, num_fulguraciones, flujo_pico_goes, rotacion_carrington,
         fecha_observacion, retrieved_at)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
''', registros)

metadatos = {
    'fuente': 'HEK (Heliophysics Event Knowledgebase) via SunPy Fido',
    'frm_regiones_activas': 'NOAA SWPC Observer',
    'evento_fulguraciones': 'FL (todos los módulos, unido por ar_noaanum)',
    'ventana_inicio': FECHA_INICIO,
    'ventana_fin': FECHA_FIN,
    'n_regiones_antes_de_cortes': str(n_antes),
    'n_regiones_despues_de_cortes': str(n_despues),
    'n_regiones_unicas': str(len(df_regiones)),
    'unidad_area': 'Mm^2 (convertido desde MSH, 1 MSH = %.4f Mm^2)' % MM2_POR_MSH,
    'unidad_flujo': 'W/m^2 (flujo de rayos X GOES en el pico de la fulguracion)',
    'generado_en': marca_tiempo,
}
cur.executemany('INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)',
                 list(metadatos.items()))

conn.commit()
conn.close()
print(f"Base de datos '{RUTA_DB}' poblada con {len(registros)} regiones activas.")
""")


# ═════════════════════════════════════════════════════════════════════════
# ETAPA 3 — Consultas SQL
# ═════════════════════════════════════════════════════════════════════════

md(r"""## Etapa 3 · Consultas SQL y agregación

Reabrimos la base (como lo haría cualquier análisis posterior, sin depender
de las variables en memoria de las etapas anteriores) y ejercitamos
`SELECT`, `WHERE` parametrizado, `GROUP BY`, `COUNT`, `AVG`, `MIN`, `MAX` y
`ORDER BY`.
""")

code(r"""conn = sqlite3.connect(RUTA_DB)

df = pd.read_sql('SELECT * FROM active_regions', conn)

# Subconjunto filtrado: regiones grandes que no sean puramente tipo Alpha
umbral_area = 100.0  # Mm^2
df_grandes_complejas = pd.read_sql(
    '''SELECT noaa_ar, area_mm2, clase_magnetica, num_fulguraciones
       FROM active_regions
       WHERE area_mm2 > ? AND clase_magnetica != ?
       ORDER BY area_mm2 DESC''',
    conn, params=(umbral_area, 'Alpha')
)

# Estadística agregada por clase magnética
df_stats = pd.read_sql(
    '''SELECT clase_magnetica,
              COUNT(*)                AS n,
              AVG(area_mm2)           AS area_media_mm2,
              MIN(area_mm2)           AS area_min_mm2,
              MAX(area_mm2)           AS area_max_mm2,
              AVG(num_fulguraciones)  AS fulguraciones_media,
              MAX(flujo_pico_goes)    AS flujo_pico_max
       FROM active_regions
       GROUP BY clase_magnetica
       ORDER BY n DESC''',
    conn
)

conn.close()

print(f"Regiones con área > {umbral_area} Mm^2 y clase != Alpha: {len(df_grandes_complejas)}")
display(df_stats)
df.describe()
""")


# ═════════════════════════════════════════════════════════════════════════
# ETAPA 4 — Dashboard estadístico
# ═════════════════════════════════════════════════════════════════════════

md(r"""## Etapa 4 · Dashboard estadístico

Construimos primero cada panel por separado (con su interpretación física
justo debajo), y al final los combinamos en un único dashboard de 3×2 que se
guarda como `dashboard.pdf`.
""")

# ── Panel 1: Pie de clase magnética ─────────────────────────────────────────
md(r"""### Panel 1 — Distribución de clases magnéticas de Hale""")
code(r"""conteo_clases = df['clase_magnetica'].value_counts()

fig, ax = plt.subplots(figsize=(6, 6))
colores = sns.color_palette('Set2', len(conteo_clases))
ax.pie(conteo_clases.values, labels=conteo_clases.index, autopct='%1.1f%%',
       colors=colores, startangle=90, pctdistance=0.8)
ax.set_title(f'Distribución de clases magnéticas de Hale (N = {len(df)} regiones)')
plt.tight_layout()
plt.show()
""")
md(r"""**Interpretación:** La mayoría de las regiones activas caen en clases
simples (`Alpha`, `Beta`), consistente con que la mayor parte del flujo
magnético emergente forma configuraciones bipolares sencillas. Las clases
complejas (`Beta-Gamma`, `Beta-Gamma-Delta`) son minoritarias pero
desproporcionadamente importantes: son las que veremos concentrar las
fulguraciones más energéticas en los paneles siguientes.
""")

# ── Panel 2: Scatter área vs fulguraciones + ajuste lineal ─────────────────
md(r"""### Panel 2 — Área vs. número de fulguraciones""")
code(r"""x = df['area_mm2'].values
y = df['num_fulguraciones'].values

pendiente, intercepto, r_valor, p_valor, error_est = stats.linregress(x, y)

fig, ax = plt.subplots(figsize=(7, 5.5))
ax.scatter(x, y, alpha=0.6, edgecolor='k', linewidth=0.3, c='steelblue')

x_ajuste = np.linspace(x.min(), x.max(), 100)
ax.plot(x_ajuste, pendiente * x_ajuste + intercepto, 'r--', lw=2,
        label=f'Ajuste lineal: R²={r_valor**2:.2f}, p={p_valor:.1e}')

ax.set_xlabel(r'Área en disco  $A$  (Mm$^2$)')
ax.set_ylabel('Número de fulguraciones asociadas')
ax.set_title('Área de la región activa vs. productividad de fulguraciones')
ax.legend()
plt.tight_layout()
plt.show()

print(f"Pendiente: {pendiente:.4f} fulguraciones/Mm^2  (R^2 = {r_valor**2:.3f}, p = {p_valor:.2e})")
""")
md(r"""**Interpretación:** Existe una correlación positiva entre el área de
la región y su número de fulguraciones: más área suele implicar más flujo
magnético libre disponible para reconectar. La dispersión es grande —el área
por sí sola no determina la actividad; la topología magnética (clase de
Hale) importa al menos tanto como el tamaño, como veremos en el panel 4.
""")

# ── Panel 3: Histograma flujo pico + ley de potencia ───────────────────────
md(r"""### Panel 3 — Distribución del flujo de rayos X pico (GOES)""")
code(r"""flujo = df.loc[df['flujo_pico_goes'] > 0, 'flujo_pico_goes'].values
log_flujo = np.log10(flujo)

fig, ax = plt.subplots(figsize=(7, 5.5))
n_bins = 12
cuentas, bordes, _ = ax.hist(log_flujo, bins=n_bins, color='darkorange',
                              edgecolor='black', alpha=0.85)
centros = 0.5 * (bordes[:-1] + bordes[1:])

# Ajuste de ley de potencia: log10(cuentas) = -alpha*log10(flujo) + b
mascara = cuentas > 0
pendiente_ley, b = np.polyfit(centros[mascara], np.log10(cuentas[mascara]), 1)

x_ajuste = np.linspace(centros.min(), centros.max(), 100)
ax.plot(x_ajuste, 10**(pendiente_ley * x_ajuste + b), 'k-', lw=2.5,
        label=fr'Ley de potencia: índice $\alpha$ = {-pendiente_ley:.2f}')

ax.set_yscale('log')
ax.set_xlabel(r'$\log_{10}$(Flujo de rayos X pico) [W/m$^2$]')
ax.set_ylabel('Número de regiones activas (escala log)')
ax.set_title('Distribución del flujo GOES pico entre regiones con fulguraciones')
ax.legend()
plt.tight_layout()
plt.show()

print(f"Índice de la ley de potencia (pendiente en log-log): {-pendiente_ley:.2f}")
""")
md(r"""**Interpretación:** La distribución de flujos pico decae como una
ley de potencia, el comportamiento estándar de la estadística de
fulguraciones solares (frecuencia decreciente con la energía liberada). El
índice obtenido es comparable en orden de magnitud a los valores
reportados en la literatura para distribuciones de flujo pico GOES
(típicamente entre ~1.5 y ~2.5, p. ej. Crosby et al. 1993), aunque nuestra
muestra —acotada a 4 meses— tiene pocas fulguraciones grandes y por tanto
el ajuste en la cola es poco confiable.
""")

# ── Panel 4: Box plot área por clase magnética ─────────────────────────────
md(r"""### Panel 4 — Área por clase magnética""")
code(r"""conteo = df['clase_magnetica'].value_counts()
clases_con_muestra = conteo[conteo >= 3].index
df_box = df[df['clase_magnetica'].isin(clases_con_muestra)]

orden_clases = (df_box.groupby('clase_magnetica')['area_mm2']
                .median().sort_values().index)

fig, ax = plt.subplots(figsize=(8, 5.5))
sns.boxplot(data=df_box, x='clase_magnetica', y='area_mm2', order=orden_clases,
            hue='clase_magnetica', palette='viridis', legend=False, ax=ax)
ax.set_xlabel('Clase magnética de Hale')
ax.set_ylabel(r'Área en disco  $A$  (Mm$^2$)')
ax.set_title('Área de la región activa por clase magnética (clases con N≥3)')
plt.xticks(rotation=20)
plt.tight_layout()
plt.show()
""")
md(r"""**Interpretación:** Las clases magnéticamente complejas (con
componentes γ o δ, que indican polaridades opuestas muy próximas) tienden a
tener áreas medianas mayores que las clases simples tipo `Alpha`. Esto es
consistente con la física de emergencia de flujo: los tubos magnéticos más
grandes y retorcidos generan configuraciones multipolares más complejas y
manchas más extensas.
""")

# ── Panel 5: Serie temporal AR por rotación de Carrington ──────────────────
md(r"""### Panel 5 — Número de regiones activas por rotación de Carrington""")
code(r"""conteo_rotacion = df.groupby('rotacion_carrington')['noaa_ar'].nunique().sort_index()

fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(conteo_rotacion.index.astype(str), conteo_rotacion.values,
       color='mediumseagreen', edgecolor='black')
ax.set_xlabel('Rotación de Carrington')
ax.set_ylabel('Número de regiones activas nuevas')
ax.set_title('Regiones activas observadas por rotación solar (Carrington)')
plt.tight_layout()
plt.show()
""")
md(r"""**Interpretación:** El número de regiones nuevas por rotación
(~27.3 días) fluctúa de una rotación a otra — refleja tanto la variabilidad
intrínseca de la emergencia de flujo cerca del máximo solar como el hecho de
que solo vemos la mitad del Sol en cada instante. No se espera una tendencia
monótona en una ventana de 4 meses; para ver el ciclo solar completo (11
años) se necesitaría extender la consulta a todo el ciclo 24.
""")

# ── Panel 6: Correlación de Pearson ────────────────────────────────────────
md(r"""### Panel 6 — Matriz de correlación de Pearson""")
code(r"""df_corr = df.copy()
df_corr['log_flujo_pico'] = np.log10(df_corr['flujo_pico_goes'].where(df_corr['flujo_pico_goes'] > 0))
df_corr['lat_abs'] = df_corr['hgs_lat_deg'].abs()

columnas_numericas = ['area_mm2', 'num_manchas', 'num_fulguraciones',
                       'log_flujo_pico', 'lat_abs']
matriz_corr = df_corr[columnas_numericas].corr(method='pearson')

fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(matriz_corr, annot=True, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1,
            square=True, cbar_kws={'label': 'Coeficiente de Pearson'}, ax=ax)
ax.set_title('Correlación entre parámetros físicos de las regiones activas')
plt.tight_layout()
plt.show()
""")
md(r"""**Interpretación:** El área y el número de manchas están fuertemente
correlacionados (esperado, ambos miden el tamaño de la región). El flujo
pico GOES correlaciona positivamente con el área y el número de
fulguraciones, reforzando que el tamaño de la región es un buen —aunque
imperfecto— predictor de la actividad energética. La correlación con la
latitud heliográfica es débil, consistente con que en esta ventana de 4
meses no cubrimos suficiente evolución del ciclo de mariposa de Spörer.
""")


# ── Dashboard combinado final ──────────────────────────────────────────────
md(r"""### Dashboard combinado

Reunimos los seis paneles anteriores en una sola figura de 3×2 y la
guardamos como `dashboard.pdf`, el entregable pedido por el enunciado.
""")

code(r"""fig, axes = plt.subplots(3, 2, figsize=(14, 16))
fig.suptitle('Dashboard: Regiones Activas Solares (HEK, 2014-01 a 2014-05)',
             fontsize=16, fontweight='bold')

# 1) Pie — clases magnéticas
axes[0, 0].pie(conteo_clases.values, labels=conteo_clases.index, autopct='%1.1f%%',
               colors=sns.color_palette('Set2', len(conteo_clases)), startangle=90)
axes[0, 0].set_title('Distribución de clases magnéticas de Hale')

# 2) Scatter — área vs fulguraciones + ajuste
axes[0, 1].scatter(x, y, alpha=0.6, edgecolor='k', linewidth=0.3, c='steelblue')
axes[0, 1].plot(x_ajuste_dummy := np.linspace(x.min(), x.max(), 100),
                 pendiente * x_ajuste_dummy + intercepto, 'r--', lw=2,
                 label=f'R²={r_valor**2:.2f}')
axes[0, 1].set_xlabel(r'Área $A$ (Mm$^2$)'); axes[0, 1].set_ylabel('N° fulguraciones')
axes[0, 1].set_title('Área vs. número de fulguraciones')
axes[0, 1].legend()

# 3) Histograma — flujo pico + ley de potencia
axes[1, 0].hist(log_flujo, bins=n_bins, color='darkorange', edgecolor='black', alpha=0.85)
axes[1, 0].plot(x_ajuste, 10**(pendiente_ley * x_ajuste + b), 'k-', lw=2.5,
                label=fr'$\alpha$={-pendiente_ley:.2f}')
axes[1, 0].set_yscale('log')
axes[1, 0].set_xlabel(r'$\log_{10}$(Flujo pico) [W/m$^2$]'); axes[1, 0].set_ylabel('N° regiones (log)')
axes[1, 0].set_title('Distribución del flujo GOES pico')
axes[1, 0].legend()

# 4) Box plot — área por clase
sns.boxplot(data=df_box, x='clase_magnetica', y='area_mm2', order=orden_clases,
            hue='clase_magnetica', palette='viridis', legend=False, ax=axes[1, 1])
axes[1, 1].set_xlabel('Clase magnética'); axes[1, 1].set_ylabel(r'Área (Mm$^2$)')
axes[1, 1].set_title('Área por clase magnética')
axes[1, 1].tick_params(axis='x', rotation=20)

# 5) Serie temporal — AR por rotación de Carrington
axes[2, 0].bar(conteo_rotacion.index.astype(str), conteo_rotacion.values,
               color='mediumseagreen', edgecolor='black')
axes[2, 0].set_xlabel('Rotación de Carrington'); axes[2, 0].set_ylabel('N° regiones nuevas')
axes[2, 0].set_title('Regiones activas por rotación solar')

# 6) Heatmap — correlación de Pearson
sns.heatmap(matriz_corr, annot=True, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1,
            square=True, cbar_kws={'label': 'Pearson r'}, ax=axes[2, 1])
axes[2, 1].set_title('Correlación entre parámetros físicos')

plt.tight_layout()
plt.savefig('dashboard.pdf', dpi=150, bbox_inches='tight')
plt.show()
print("Dashboard guardado en 'dashboard.pdf'")
""")


# ═════════════════════════════════════════════════════════════════════════
# ETAPA 5 — Interpretación y resumen
# ═════════════════════════════════════════════════════════════════════════

md(r"""## Etapa 5 · Resumen de hallazgos

**Hallazgos principales:**

1. Las clases magnéticas simples (`Alpha`, `Beta`) dominan numéricamente el
   catálogo, pero las clases complejas (`Beta-Gamma-Delta`) concentran una
   fracción desproporcionada del área y del flujo de rayos X pico.
2. El área en disco correlaciona positivamente con el número de
   fulguraciones asociadas (R² moderado) — el tamaño ayuda a predecir la
   actividad, pero no la determina por completo.
3. La distribución del flujo GOES pico sigue aproximadamente una ley de
   potencia, consistente con la fenomenología estándar de fulguraciones
   solares.
4. La complejidad magnética (clase de Hale) se asocia con áreas medianas
   mayores, apoyando el vínculo físico entre topología magnética y tamaño
   de la región.
5. En una ventana de solo 4 meses no se observa una tendencia clara del
   número de regiones nuevas por rotación de Carrington ni de la latitud
   heliográfica — se necesitaría el ciclo solar completo para verlo.

**Limitaciones:**

- El catálogo NOAA/SWPC solo ve el hemisferio visible desde la Tierra:
  regiones que emergen y decaen en la cara oculta nunca se registran.
- La ventana de 4 meses es corta para las conclusiones sobre evolución del
  ciclo solar (ítems 5 y el panel de correlación con latitud).
- El ajuste de ley de potencia en el panel 3 usa pocos bins (muestra
  limitada) y no aplica un método riguroso de máxima verosimilitud —una
  extensión natural sería usar la librería `powerlaw` de Alstott et al.
  (2014) para un ajuste estadísticamente más robusto.

**Posibles análisis de seguimiento:**

- Extender la ventana de consulta al ciclo solar 24 completo (2008–2019)
  para estudiar el diagrama de mariposa y la relación entre fase del ciclo
  y complejidad magnética.
- Incorporar magnetogramas de SHARP (HMI) para calcular directamente el
  flujo magnético no potencial en vez de usar solo la clase de Hale como
  proxy categórico.
- Modelar la probabilidad de fulguración como función logística de área y
  clase magnética (similar a los esquemas operacionales de predicción de
  clima espacial, p. ej. McIntosh o el índice de Falconer).

**Referencias:**

- Hoeksema, J. T. et al. (2014), *The Helioseismic and Magnetic Imager
  (HMI) Vector Magnetic Field Pipeline*, Solar Physics.
- Crosby, N. B., Aschwanden, M. J., & Dennis, B. R. (1993), *Frequency
  distributions and correlations of solar X-ray flare parameters*, Solar
  Physics 143, 275.
- The SunPy Community (2020), *The SunPy Project: Open Source Development
  and Status of the Version 1.0 Core Package*, ApJ 890, 68.
- NOAA Space Weather Prediction Center — documentación del esquema de
  clasificación magnética de Hale y del catálogo de regiones activas.
- Documentación del HEK: https://www.lmsal.com/hek/
""")


# ═════════════════════════════════════════════════════════════════════════
# Escribir el cuaderno
# ═════════════════════════════════════════════════════════════════════════

cuaderno['cells'] = celdas

with open('solar_active_regions.ipynb', 'w') as archivo:
    nbf.write(cuaderno, archivo)

print("✅  solar_active_regions.ipynb escrito exitosamente!")
