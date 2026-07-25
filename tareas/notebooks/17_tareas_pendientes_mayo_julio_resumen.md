# Tareas pendientes mayo-julio 2026

Este resumen fue generado por `notebooks/17_tareas_pendientes_mayo_julio.py`.

## Optimizacion
- Temperatura optima: T = 30793.74 K, P = 3.5140e+08 W m^-2.
- Disco de acrecion: todos los inicios convergen al maximo esperado r/r0 ~= 2, theta ~= pi/3.

## PDE elipticas
- Gauss-Seidel: 163 iteraciones, error aprox. 0.0987%.
- SOR lambda=1.5: 79 iteraciones, error aprox. 0.0964%.
- Potencial minimo SOR: -8.3801e+06 m^2 s^-2.

## FEM 1D
- n=10: error relativo maximo 1.344e-02, |g(R)| ~= 8.475e+03 m s^-2.
- n=20: error relativo maximo 4.194e-03, |g(R)| ~= 8.184e+03 m s^-2.
- n=40: error relativo maximo 1.255e-03, |g(R)| ~= 8.015e+03 m s^-2.
- Valor teorico del modelo: |g(R)| = 7.828e+03 m s^-2; valor solar de referencia = 274.0 m s^-2.

## MHD
- Tiempo para recorrer 1 AU a 400 km/s: 4.329 dias.
- Velocidad de Alfven del ejemplo: 892.062 km/s.

## FITS
- Fondo mediano removido: 9658.000; sigma robusta: 3881.447.
- FITS modificado: `notebooks/17_tareas_pendientes_outputs/21_horsehead_modificado.fits`.

## SQL
- Shows con rating 10.0: 27.
- Episodios de Black Mirror: 22.
- Shows Sci-Fi: 986.
- Mejor Horror: O Kapoios (10.0).
- Shows Animation: 3104.

## Machine Learning
- Hubble: pendiente=3.5500e-04, R2=0.767.
- SMBH M-sigma: pendiente=2.925, R2=0.619.
- Bolshoi: R2 test=0.821, RMSE test=0.133.
- DBSCAN galaxias/fotometria: 3 clusters y 5302 puntos de ruido sobre 20000.

## Figuras
- `notebooks/17_tareas_pendientes_figuras/17_temperatura_optima.png`
- `notebooks/17_tareas_pendientes_figuras/17_disco_acrecion_gradiente.png`
- `notebooks/17_tareas_pendientes_figuras/18_potencial_disco_protoplanetario.png`
- `notebooks/17_tareas_pendientes_figuras/18_fuerza_radial_disco.png`
- `notebooks/17_tareas_pendientes_figuras/19_fem_poisson_gravitacional_1d.png`
- `notebooks/17_tareas_pendientes_figuras/20_mhd_adveccion_1d.png`
- `notebooks/17_tareas_pendientes_figuras/20_onda_alfven.png`
- `notebooks/17_tareas_pendientes_figuras/21_fits_horsehead_modificado.png`
- `notebooks/17_tareas_pendientes_figuras/22_ml_regresiones.png`
- `notebooks/17_tareas_pendientes_figuras/23_clustering_galaxias_dbscan.png`
