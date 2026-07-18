---
title: "Regiones Activas Solares: De la Consulta al Dashboard"
subtitle: "Proyecto 3 — Astrofísica Computacional · Opción A (HEK / SunPy)"
author: "Carlos Alberto Martínez Sibaja"
date: "Julio de 2026"
geometry: margin=2.5cm
fontsize: 11pt
---

## 1. Introducción

Una región activa (AR, *Active Region*) es una zona de la fotosfera solar
donde un tubo de flujo magnético emerge desde la zona de convección y forma
manchas solares. NOAA clasifica su topología magnética con el esquema de
Hale/Mt. Wilson (α, β, γ, δ y combinaciones), y esta clasificación es uno de
los mejores predictores operacionales de si una región producirá una
fulguración solar.

Este proyecto construye un pipeline completo de datos astronómicos —consulta
real, base de datos relacional, análisis SQL y dashboard estadístico— usando
el catálogo de regiones activas y fulguraciones del *Heliophysics Event
Knowledgebase* (HEK), consultado con SunPy. Se eligió esta opción (Opción A
del enunciado) por continuidad temática con los proyectos anteriores del
curso, ambos de física solar (viento solar de Parker; reconexión magnética
Hall MHD).

**Pregunta central:** ¿qué tan bien predicen el tamaño y la clase magnética
de una región activa su productividad de fulguraciones y la energía de
rayos X liberada?

## 2. Métodos

### 2.1 Fuente de datos y consulta

Se consultó el HEK vía `sunpy.net.hek`/`Fido` para dos tipos de evento en la
ventana 2014-01-01 a 2014-05-01 (cuatro meses cerca del máximo del ciclo
solar 24):

- **`AR`** (regiones activas), filtrando `frm_name == 'NOAA SWPC Observer'`
  — el único módulo del HEK que reporta número NOAA y clase magnética.
- **`FL`** (fulguraciones), sin filtro de módulo en la consulta inicial,
  para poder comparar la cobertura de los distintos catálogos automáticos.

La consulta devolvió 890 filas diarias de regiones activas y 10859 filas de
fulguraciones. Por el costo de red (varios minutos por mes consultado), el
resultado crudo se cachea en `data/*.csv`; el cuaderno vuelve a consultar el
HEK en vivo si el caché no existe.

### 2.2 Limpieza de datos: tres hallazgos no triviales

Inspeccionar los datos crudos (no solo la documentación del HEK) reveló tres
problemas que, sin corregir, habrían arruinado silenciosamente el análisis:

1. **Área en unidades equivocadas.** `area_atdiskcenter` no está en
   "millonésimas de hemisferio solar" (MSH, la unidad que NOAA reporta
   originalmente) sino ya en km² —un parámetro derivado que calcula el
   propio HEK—. Se detectó porque los valores estaban cuantizados en pasos
   de ≈3.04×10⁷, exactamente 10 MSH expresadas en km². Aplicar el factor de
   conversión MSH→Mm² sobre un valor que ya estaba en km² habría inflado
   las áreas por un factor ~3×10⁶.
2. **Clase magnética vacía en el campo "canónico".** `ar_noaaclass` (el
   campo que en teoría trae la clase de Hale) está vacío en el 100% de las
   filas de este catálogo. La clasificación real vive en `ar_mtwilsoncls`,
   pero con un defecto de origen: el guion entre componentes a veces se
   pierde (`"BETAAGAMMA-DELTA"` en vez de `"BETA-GAMMA-DELTA"`). Se
   reconstruyó la clase canónica con una expresión regular que extrae los
   tokens válidos (Alpha/Beta/Gamma/Delta) y colapsa repeticiones
   consecutivas.
3. **El flujo y el número NOAA nunca coexisten.** `ar_noaanum` (número de
   región) y `fl_peakflux` (flujo numérico, calculado por el módulo
   automático *Flare Detective*) nunca aparecen en la misma fila —son
   catálogos independientes del mismo evento físico—. Las filas del módulo
   `SWPC` sí traen ambos, pero el flujo solo como texto GOES (`fl_goescls`,
   p. ej. `"M9.9"`), que se convirtió a W/m² con la escala estándar
   ($10^{-8}$ a $10^{-4}$ W/m² para clases A–X).

Ninguno de estos tres problemas producía un error de Python: los tres se
detectaron verificando que los resultados intermedios fueran físicamente
razonables, no solo que el código corriera sin excepciones.

### 2.3 De observaciones diarias a un catálogo por región

NOAA reporta cada AR una vez por día mientras cruza el disco visible. Se
agregó a una fila por región NOAA tomando: el área pico observada, las
coordenadas heliográficas de Stonyhurst en el instante más cercano al
meridiano central, la clase magnética modal, el número de manchas pico, y
—uniendo con la tabla de fulguraciones por número NOAA— el conteo de
fulguraciones y el flujo GOES pico asociado.

### 2.4 Base de datos y consultas SQL

Se diseñó una base SQLite (`astro_project.db`) con una tabla
`active_regions` (una fila por región, 115 en total) y una tabla
`metadata` autodescriptiva (fuente, ventana de tiempo, cortes de calidad
aplicados, convenciones de unidades). Se ejercitaron consultas con
`SELECT`, `WHERE` parametrizado, `GROUP BY`, `COUNT`, `AVG`, `MIN`, `MAX` y
`ORDER BY` para obtener estadísticas agregadas por clase magnética.

## 3. Resultados

De las 890 filas diarias, ninguna se descartó en los cortes de calidad
(número NOAA, clase magnética y coordenadas válidas), resultando en **115
regiones activas distintas**, 71 de ellas con al menos una fulguración
asociada.

![Dashboard de 6 paneles: distribución de clases magnéticas, área vs.
fulguraciones, distribución del flujo GOES pico, área por clase magnética,
regiones por rotación de Carrington, y correlación de
Pearson.](figures/dashboard-1.png){width=95%}

**Distribución de clases magnéticas (Panel 1).** En esta muestra predominan
las configuraciones magnéticas **complejas**: Beta-Gamma-Delta (43.5%) y
Beta-Gamma (30.4%) suman en conjunto cerca del 74% de las regiones, mientras
que las clases simples (Alpha, Beta) representan apenas ~7.8%. Este
resultado no debe generalizarse a toda la población solar: la ventana
corresponde a cuatro meses cercanos al máximo del ciclo solar 24, y se usó
la clase *dominante* durante el tránsito visible de cada región (no cada
observación diaria individual).

**Área vs. fulguraciones (Panel 2).** Existe una asociación positiva y
estadísticamente significativa entre el área de la región y su número de
fulguraciones (R² = 0.71, p = 2.8×10⁻³²). Esta relación es **descriptiva**,
no un modelo predictivo ni una prueba de causalidad: el número de
fulguraciones es una variable de conteo con pocos puntos extremos que pesan
mucho en el ajuste, y no se controló por clase magnética ni por otros
factores. La dispersión es grande: el área por sí sola no determina la
actividad.

**Flujo GOES pico (Panel 3).** La distribución presenta una disminución
aproximadamente compatible con una tendencia de ley de potencia, pero el
ajuste (α ≈ 0.33) es **exploratorio** y no equivale a estudiar la
distribución completa de fulguraciones: con solo 71 regiones, 12 intervalos
y el flujo *máximo por región* (una estadística de valores extremos, no la
distribución completa de eventos), el índice obtenido **no** debe
compararse directamente con los reportados en la literatura para la
distribución completa (Crosby et al. 1993, ~1.5–2.5).

**Área por clase magnética (Panel 4).** Las clases con componente δ
(Beta-Gamma-Delta) muestran áreas medianas mayores y una cola hacia
valores extremos (hasta ~4800 Mm²), consistente con que las
configuraciones magnéticas más complejas provienen de tubos de flujo más
grandes y retorcidos.

**Regiones por rotación de Carrington (Panel 5).** El número de regiones
nuevas por rotación (~27.3 días) varía entre 19 y 27 en la ventana
estudiada. Las rotaciones de los extremos (2145 y 2149, en gris en la
figura) están cubiertas solo **parcialmente**, porque la consulta empieza
el 2014-01-01 y termina el 2014-05-01 cortando a la mitad esas dos
rotaciones — parte de la variación observada es este efecto de borde del
muestreo, no necesariamente una señal física del ciclo solar.

**Correlación de Pearson (Panel 6).** El área y el número de manchas están
fuertemente correlacionados (r = 0.93, esperado). El flujo pico GOES
correlaciona moderadamente con el área (r = 0.54) y fuertemente con el
número de fulguraciones (r = 0.62). La correlación con la latitud
heliográfica es débil (|r| < 0.15), consistente con no cubrir suficiente
evolución del ciclo de mariposa de Spörer en 4 meses.

**Validación con un evento histórico conocido:** la región con el flujo
GOES pico más alto del catálogo es NOAA AR 11990 (4.9×10⁻⁴ W/m², clase
X4.9), que coincide exactamente con la fulguración X4.9 real ocurrida el
25 de febrero de 2014 — una verificación cruzada útil de que el pipeline
recuperó datos físicamente correctos.

## 4. Conclusiones

1. El tamaño de una región activa se asocia positiva y significativamente
   con su productividad de fulguraciones (relación descriptiva, no causal
   ni predictiva), pero con dispersión grande: el tamaño ayuda a explicar
   la actividad, no la determina.
2. La complejidad magnética (clases con componente δ) se asocia con áreas
   medianas mayores y con los flujos GOES más energéticos del catálogo,
   apoyando el vínculo físico bien establecido entre topología magnética y
   liberación de energía.
3. Los datos "crudos" de un archivo público como el HEK requieren
   verificación activa: tres problemas de unidades, columnas vacías y
   granularidad de codificación se encontraron inspeccionando los valores,
   no la documentación, y ninguno habría producido un error de programa.
4. Una ventana de 4 meses es adecuada para relaciones estructurales
   (área–clase–actividad) pero insuficiente para tendencias de ciclo solar
   (variación con rotación de Carrington o latitud), y además introduce un
   efecto de borde: las rotaciones de Carrington en los extremos de la
   ventana quedan parcialmente muestreadas.

**Limitaciones y trabajo futuro:** extender la consulta al ciclo solar 24
completo permitiría estudiar el diagrama de mariposa; usar directamente la
distribución completa de fulguraciones (no el máximo por región) daría un
ajuste de ley de potencia comparable con la literatura; e incorporar
magnetogramas SHARP (HMI) permitiría reemplazar la clase de Hale
categórica por una medida continua de flujo magnético no potencial.

## 5. Referencias

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
