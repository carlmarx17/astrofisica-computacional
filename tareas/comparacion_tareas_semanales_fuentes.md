# Comparacion de tareas semanales con fuentes oficiales

Fecha de revision: 2026-07-25

Este cruce usa la lista con codigos `S...` enviada por el usuario y las fuentes del repositorio
`Hypnus1803/ComputationalAstrophysics_2026-I`. Se revisaron los notebooks locales en `notebooks/`
y, cuando hizo falta, se descargaron temporalmente PDFs/notebooks fuente en `/tmp/tareas_fuentes`.

## Resumen ejecutivo

- Encontradas en notebooks locales, incluyendo entregas consolidadas: 19.
- Parciales o agrupadas, pero no exactamente completas: 3.
- No encontrada como entrega resuelta: 1.

La tarea claramente faltante es:

- `S180402` - Pulsacion radial de una estrella.

Las tareas parciales que conviene completar o separar mejor son:

- `S300501` - Busquedas en SDSS.
- `S300502` - Los cuatro ejercicios de Astroquery2.
- `S060601` - Los otros tres ejercicios de regresion lineal con ML.

## Tabla por codigo

| Codigo | Tarea | Estado local | Notebook/evidencia local | Fuente revisada |
|---|---|---:|---|---|
| S280206 | Regla Simpson 3/8 funcion propia | Encontrada | `notebooks/08_integracion_trapezoide_y_simpson.ipynb` | `Integracion.pdf` |
| S140301 | Biseccion pelota flotante | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` | `Roots_searching.pdf` |
| S140302 | Newton-Raphson pelota flotante | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` | `Roots_searching.pdf` |
| S140303 | Secante pelota flotante | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` | `Roots_searching.pdf` |
| S140304 | Newton-Raphson ecuaciones simultaneas no lineales | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` | `Roots_searching.pdf` |
| S280301 | Biseccion para decaimiento | Encontrada | `notebooks/10_euler_ode_y_biseccion_decaimiento.ipynb` | `Roots_searching.pdf` |
| S280302 | ODE usando metodo de Euler | Encontrada | `notebooks/10_euler_ode_y_biseccion_decaimiento.ipynb` | `odes_and_euler.pdf` |
| S280401 | Tareas Runge Kutta | Encontrada | `notebooks/11_runge_kutta_tareas.ipynb` | `runge_kutta.pdf` |
| S180401 | Perfil de temperatura zona radiativa solar | Encontrada | `notebooks/13_diferencias_finitas_temperatura_solar.ipynb` | `odes_2_raw.pdf` |
| S180402 | Pulsacion radial de una estrella | No encontrada | No hay notebook local con la aplicacion astronomica de pulsacion radial. | `odes_2_raw.pdf` |
| S250401 | Temperatura optima de emision estelar | Encontrada | `notebooks/15_seccion_aurea_temperatura_optima_emision_estelar.ipynb`; tambien aparece en `17`. | `Optimizacion.pdf` |
| S250402 | Maximizar la luminosidad de un disco de acrecion | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` | `Optimizacion.pdf` |
| S250403 | Problemas anteriores usando metodo de la gradiente | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` | `Optimizacion.pdf` |
| S090501 | Enfriamiento termico corteza estrella de neutrones | Encontrada | `notebooks/16_diferencias_finitas_enfriamiento_corteza_estrella_neutrones.ipynb` | `derivadas_parciales01.pdf` |
| S160501 | Potencial gravitacional en disco protoplanetario | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` | `derivadas_parciales_ellipticas.pdf` |
| S160502 | Adveccion 1D | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` | `MHD.pdf` |
| S160502 | FEM para Poisson gravitacional 1D | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` | `derivadas_parciales_ellipticas.pdf` |
| S230501 | Recreacion de ondas Alfven | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` | `MHD.pdf` |
| S300501 | Modificar los FITS del notebook | Encontrada, consolidada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` y `notebooks/17_tareas_pendientes_outputs/21_horsehead_modificado.fits` | `01. FITSImage01.ipynb` |
| S300501 | Busquedas en SDSS | Parcial | `17` menciona estructura SDSS/Astroquery sin depender de red, pero no contiene busquedas SDSS ejecutadas ni resultados verificables. | `01. SQLExample.ipynb` |
| S300502 | Los cuatro ejercicios de notebook | Parcial | `17` cubre FITS/SQL/SDSS de forma consolidada, pero no resuelve explicitamente los cuatro ejercicios de `Astroquery2`. | `04. Astroquery2.ipynb` |
| S060601 | Los otros tres ejercicios de regresion lineal con ML | Parcial | `17` contiene Hubble, SMBH y Bolshoi; la fuente oficial separa los otros tres como `02`, `03`, `04` de ScikitRegression. Falta evidencia de Ridge/LASSO y del notebook 04 completo. | `ML I. Regression/02-04` |
| S270602 | Clustering de galaxias con datos fotometricos reales | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` y figura `23_clustering_galaxias_dbscan.png` | `Clustering_AstroML_DatosReales.ipynb` |

## Detalles fuente utiles para completar faltantes

### S180402 - Pulsacion radial de una estrella

La fuente `odes_2_raw.pdf` pide:

- Transformar una ODE de segundo orden a dos ODEs de primer orden.
- Resolver con RK4 con `h = 0.01`.
- Graficar `x(t)` y `v(t)` para observar oscilaciones amortiguadas.

Modelo indicado por la fuente:

```text
x(t) = R(t) - R0
d2x/dt2 + 2 xi omega0 dx/dt + omega0^2 x = 0
omega0 = 2 pi / P
P = 5 dias ~= 4.32e5 s
xi = 0.1
R0 ~= 42.7 R_sun
```

Estado local: falta una entrega explicita. El notebook `14_metodo_disparo_y_odes_orden_superior.ipynb`
cubre metodo de disparo y ODEs de orden superior, pero no esta aplicacion astronomica.

### S300502 - Cuatro ejercicios de Astroquery2

La fuente `04. Astroquery2.ipynb` pide:

1. Identificar al menos tres lineas del espectro de NGC5406 y calcular su redshift.
2. Descargar imagenes FITS en filtros `u,g,r,z,i`, graficar logaritmo del flujo en cinco paneles y explicar diferencias.
3. Calcular el perfil radial de flujo de NGC5406 con 10 lineas desde el centro y proponer una funcion de ajuste.
4. Repetir el procedimiento para `SDSS J013755.71+010004.9` y explicar por que se ve diferente de NGC5406.

Estado local: `17` no contiene esos cuatro ejercicios como resultados separados.

### S060601 - Otros tres ejercicios de regresion lineal con ML

La fuente oficial contiene cuatro notebooks:

- `01. ScikitRegression01.ipynb`: ajuste lineal de datos de Hubble.
- `02. ScikitRegression02.ipynb`: relacion `M_BH - sigma_*` y `M_BH - FWHM` para agujeros negros supermasivos.
- `03. ScikitRegression03.ipynb`: ajuste multilineal, Ridge y LASSO con los datos SMBH.
- `04. ScikitRegression04.ipynb`: ajuste lineal con `data_sample.csv`.

Estado local: `17` cubre Hubble, SMBH y Bolshoi, pero no deja evidencia completa de los notebooks `02`, `03` y `04`
segun la estructura oficial.

## Correcciones de datos

- `S180402` aparece con fecha final `2/4/2026`, anterior al inicio `18/4/2026`. Probablemente debe ser `25/4/2026` o `2/5/2026`.
- El codigo `S160502` esta duplicado para `Adveccion 1D` y `FEM para Poisson gravitacional (1D)`. Conviene renombrar uno si se usa como identificador unico.
- El codigo `S300501` esta duplicado para `Modificar FITS` y `Busquedas en SDSS`. Conviene renombrar uno si se usa como identificador unico.
