# 🔆 Regiones Activas Solares — De la Consulta al Dashboard

> *"¿Por qué unas manchas solares producen fulguraciones X y otras nunca
> hacen nada? La respuesta está en su topología magnética, no solo en su
> tamaño. Aquí lo comprobamos con datos reales del HEK."*

---

## ¿Qué es este proyecto?

Este proyecto construye un **pipeline completo de datos astronómicos**:
consulta real al *Heliophysics Event Knowledgebase* (HEK) vía SunPy,
ingesta en una base de datos SQLite normalizada, análisis estadístico con
SQL + pandas, y un dashboard de 6 paneles con interpretación física.

Elegimos la **Opción A** del enunciado (Regiones Activas Solares) porque
continúa la línea de física solar de los proyectos anteriores del curso
(viento solar de Parker, reconexión magnética Hall MHD).

---

## La Física en Dos Oraciones

Una región activa es donde un tubo de flujo magnético emerge desde la zona
de convección y forma manchas solares. NOAA clasifica su topología con el
**esquema de Hale** (α, β, γ, δ y combinaciones como β-γ-δ); cuanto más
compleja la clase —más polaridades opuestas entrelazadas—, más energía libre
hay disponible para liberarse en una fulguración.

---

## Estructura de Archivos

```
proyecto3/
│
├── project3.pdf                    ← Enunciado del proyecto
├── generar_notebook.py             ← Construye el cuaderno programáticamente
├── solar_active_regions.ipynb      ← El reporte principal (generado)
├── astro_project.db                ← Base SQLite poblada por el cuaderno
├── dashboard.pdf                   ← Dashboard de 6 paneles (generado)
│
├── data/
│   ├── fetch_raw_data.py           ← Descarga cruda del HEK (AR + FL, 4 meses)
│   ├── ar_raw.csv                  ← Caché: regiones activas (NOAA SWPC Observer)
│   └── fl_raw.csv                  ← Caché: fulguraciones (todas, con ar_noaanum)
│
├── report/
│   ├── report.md                   ← Fuente del reporte (introducción, métodos, resultados, conclusiones)
│   ├── report.pdf                  ← Reporte corto (4 páginas), generado con pandoc + weasyprint
│   ├── style.css                   ← Estilo de márgenes/tipografía para el PDF
│   └── figures/dashboard-1.png     ← Figura del dashboard embebida en el reporte
│
├── requirements.txt                ← Dependencias exactas (pip freeze verificado)
└── README.md                       ← Estás aquí
```

---

## Fuente de Datos

| | |
|---|---|
| **Servicio** | HEK (Heliophysics Event Knowledgebase), vía `sunpy.net.hek` / `Fido` |
| **Eventos consultados** | `AR` (regiones activas, filtro `frm_name == 'NOAA SWPC Observer'`) y `FL` (fulguraciones, todos los módulos) |
| **Ventana de tiempo** | 2014-01-01 a 2014-05-01 (4 meses, cerca del máximo del ciclo solar 24) |
| **Unión** | Las fulguraciones se asocian a su región activa por `ar_noaanum` |

**¿Por qué filtrar `frm_name`?** El HEK cataloga el mismo tipo de evento
(`AR`) con varios módulos automáticos de reconocimiento (SPoCA, HMI SHARP,
NOAA SWPC Observer). Solo el reporte oficial de **NOAA SWPC** trae el
número de región NOAA y la clase magnética de Hale que pide el enunciado;
los otros dos no tienen esas columnas.

**¿Por qué cachear en `data/*.csv`?** Cada mes de consulta al HEK tarda
minutos (cuello de botella del servidor, no del código). El cuaderno
detecta el caché y lo usa si existe; si lo borras, vuelve a consultar el
HEK en vivo con el mismo código de `data/fetch_raw_data.py` reproducido en
la Etapa 1 del cuaderno.

---

## Diseño de la Base de Datos

Una fila por región activa NOAA (no por observación diaria — ver
"Agregación" más abajo):

```sql
CREATE TABLE active_regions (
    id                   INTEGER PRIMARY KEY AUTOINCREMENT,
    noaa_ar              INTEGER NOT NULL UNIQUE,
    hgs_lon_deg          REAL,     -- Longitud heliográfica de Stonyhurst [°]
    hgs_lat_deg          REAL,     -- Latitud heliográfica de Stonyhurst [°]
    area_mm2             REAL,     -- Área pico en disco [Mm²]
    clase_magnetica      TEXT,     -- Clase de Hale dominante (α/β/γ/δ)
    num_manchas          INTEGER,  -- Número de manchas pico
    num_fulguraciones    INTEGER,  -- Fulguraciones asociadas (0 si ninguna)
    flujo_pico_goes      REAL,     -- Flujo de rayos X pico [W/m²] (NULL si no hubo)
    rotacion_carrington  INTEGER,  -- Rotación de Carrington de 1ª observación
    fecha_observacion    TEXT,
    retrieved_at         TEXT
);

CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT);  -- auto-documentación
```

### De observaciones diarias a un catálogo por región

NOAA reporta cada región activa **una vez por día** mientras cruza el disco
visible, así que el HEK devuelve muchas filas por región. Agregamos así:

| Parámetro | Cómo se calcula |
|---|---|
| Área [Mm²] | Máximo observado durante el tránsito (convertido de MSH a Mm²) |
| Coordenadas heliográficas | Las del instante más cercano al meridiano central (`\|hgs_x\|` mínimo) |
| Clase magnética | Moda (clasificación más frecuente) durante el tránsito |
| N° de manchas | Máximo observado |
| N° fulguraciones / flujo pico | Unión con la tabla `FL` por `ar_noaanum`; 0 fulguraciones es información física real, no un dato faltante |

**Conversión de unidades:** NOAA/SWPC reporta el área en "millonésimas de
hemisferio solar" (MSH). La conversión usada es:

$$1\ \text{MSH} = \frac{2\pi R_\odot^2}{10^6} \approx 3.04\ \text{Mm}^2$$

---

## Contenido del Cuaderno (`solar_active_regions.ipynb`)

| Etapa | Qué cubre |
|---|---|
| **1 — Adquisición** | Consulta HEK (AR + FL), cortes de calidad, discusión de sesgos de selección |
| **2 — Base de datos** | Esquema SQLite, agregación por región, tabla `metadata` auto-documentada |
| **3 — Consultas SQL** | `SELECT`, `WHERE` parametrizado, `GROUP BY`, `COUNT/AVG/MIN/MAX` |
| **4 — Dashboard** | 6 paneles: pie de clases, área vs. fulguraciones (ajuste lineal), histograma de flujo pico (ley de potencia), boxplot de área por clase, serie temporal por rotación de Carrington, heatmap de correlación de Pearson |
| **5 — Interpretación** | Hallazgos, limitaciones, análisis de seguimiento, referencias |

---

## Cómo Correr el Proyecto

```bash
# 1. Crear el entorno virtual e instalar dependencias exactas
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# 2. (Opcional) Regenerar el cuaderno desde el código fuente
python generar_notebook.py

# 3. Ejecutarlo de punta a punta
jupyter nbconvert --to notebook --execute --inplace solar_active_regions.ipynb

# 4. O abrirlo interactivamente
jupyter notebook solar_active_regions.ipynb
```

Si `data/ar_raw.csv` y `data/fl_raw.csv` ya existen, el paso 3 tarda
segundos. Si los borras, la Etapa 1 vuelve a consultar el HEK en vivo
(~15 minutos para la ventana de 4 meses).

**Para regenerar `report/report.pdf`** (no requiere las dependencias de
`requirements.txt`, solo `pandoc` y `weasyprint`):

```bash
pip install --user weasyprint
cd report && pandoc report.md -o report.pdf --pdf-engine=weasyprint -c style.css
```

---

## Resultados Clave de un Vistazo

Corrida real con datos del HEK (2014-01-01 a 2014-05-01):

| Cantidad | Valor |
|---|---|
| Filas diarias recuperadas (AR) | 890 (0 descartadas en cortes de calidad) |
| Regiones activas NOAA distintas | 115 |
| Regiones con ≥1 fulguración asociada | 71 |
| Área en disco: mínima / media / máxima | 30.4 / 492.4 / 4809.0 Mm² |
| Clases magnéticas complejas (Beta-Gamma + Beta-Gamma-Delta) | ~74% de la muestra (Alpha+Beta: ~7.8%) — no generalizable, ver nota abajo |
| Asociación área–fulguraciones (descriptiva) | R² = 0.71 (p = 2.8×10⁻³²) |
| Flujo GOES pico más alto del período | 4.9×10⁻⁴ W/m² (X4.9, NOAA AR 11990 — coincide con la fulguración X4.9 real del 25 de febrero de 2014) |
| Correlación área–número de manchas | Pearson r = 0.93 |

**Notas de interpretación importantes** (ver el cuaderno para el detalle):
el predominio de clases complejas es específico de esta ventana de 4 meses
cerca del máximo del ciclo solar 24 y no debe generalizarse; la regresión
área–fulguraciones es descriptiva, no un modelo causal ni predictivo; el
ajuste de ley de potencia del flujo GOES (Panel 3) es exploratorio (usa
solo el flujo máximo por región, no la distribución completa de
fulguraciones); y las rotaciones de Carrington en los bordes de la consulta
(2145, 2149) están parcialmente muestreadas, por lo que su menor conteo
puede ser un efecto de borde y no una señal física real.

**Tres hallazgos de calidad de datos** encontrados inspeccionando los CSV
crudos (documentados con detalle en el cuaderno):

1. `area_atdiskcenter` del HEK ya viene en km² (no en "millonésimas de
   hemisferio solar" como reporta NOAA originalmente) — se detectó por la
   cuantización de los valores en pasos de exactamente 10 MSH.
2. El campo "canónico" de clase magnética del HEK (`ar_noaaclass`) está
   completamente vacío; hubo que usar `ar_mtwilsoncls` con una limpieza por
   expresión regular (guiones faltantes entre componentes).
3. `ar_noaanum` y el flujo numérico de fulguraciones nunca coexisten en la
   misma fila del HEK; el flujo se reconstruyó desde la clase GOES en texto
   (`fl_goescls`) con la escala estándar A–X.

---

*Repositorio del curso Astrofísica Computacional · 2026-I · OAN*
