# Comparacion de tareas contra notebooks

Fecha de revision: 2026-07-25

Nota: tambien existe un cruce mas preciso por codigo semanal y fuente oficial en
`comparacion_tareas_semanales_fuentes.md`.

## Resumen

- Notebooks encontrados: `notebooks/01_...ipynb` a `notebooks/17_...ipynb`.
- Indice existente revisado: `tareas_resueltas.md`.
- Muchas tareas diarias no estan como archivo separado, pero si aparecen dentro de notebooks semanales.
- Las tareas mas claramente cubiertas estan entre febrero y mayo, mas el consolidado `17_tareas_pendientes_mayo_julio.ipynb`.
- Los huecos mas importantes estan en regresion exponencial, algunos ejercicios de optimizacion, orbitas, cilindro de pared gruesa y clasificacion ML 2.

## Leyenda

- `Encontrada`: aparece de forma explicita por nombre, titulo o implementacion.
- `Parcial/probable`: hay contenido relacionado, pero no coincide completamente con el nombre de la lista.
- `No encontrada`: no encontre una coincidencia clara en notebooks, scripts, resumenes o figuras.
- `Fuera de notebook`: evidencia en el repo, pero no como notebook.

## Comparacion por tarea

| Fecha | Tarea de la lista | Estado | Evidencia |
|---|---|---:|---|
| 21/2-28/2 | Crear repositorio | Fuera de notebook | El repo existe y contiene `README.md`, `tareas_resueltas.md` y notebooks. |
| 21/2-28/2 | Galaxia espiral | Encontrada | `notebooks/01_galaxia_espiral_y_04_ejercicios.ipynb` |
| 21/2-28/2 | Ejercicio de la clase del 14 de febrero | Encontrada | `notebooks/02_clase_14_febrero_y_repaso.ipynb` |
| 21/2-28/2 | Tarea 14 de febrero | Encontrada | `notebooks/03_tarea_14_febrero_03_ejercicios.ipynb` |
| 21/2-28/2 | More differentation examples | Parcial/probable | `notebooks/04_picos_verdes_en_grafica.ipynb` y ejercicios de diferencias finitas en `03`. |
| 21/2-28/2 | Exercise01 | Encontrada | `notebooks/05_todo_exercises01.ipynb` |
| 21/2-28/2 | Exercise python / promedio de temperatura | Parcial/probable | `notebooks/05_todo_exercises01.ipynb` contiene datos de temperatura, media movil y clases; `notebooks/06_respuestas_05_preguntas_exercises_python.ipynb` contiene las 5 preguntas de Exercises Python. |
| 28/2 | Integracion del cohete | Encontrada | `notebooks/07_integracion_cohete_y_riemann.ipynb` y `notebooks/08_integracion_trapezoide_y_simpson.ipynb` |
| 28/2 | Metodo del trapecio simple | Encontrada | `notebooks/08_integracion_trapezoide_y_simpson.ipynb` |
| 28/2 | Metodo del trapecio compuesto | Encontrada | `notebooks/08_integracion_trapezoide_y_simpson.ipynb` |
| 28/2 | Metodo de Simpson | Encontrada | `notebooks/08_integracion_trapezoide_y_simpson.ipynb` |
| 28/2-7/3 | Integracion analitica del cohete | Encontrada | `notebooks/07_integracion_cohete_y_riemann.ipynb` |
| 28/2-7/3 | Integrar usando Riemann y=x^3 | Encontrada | `notebooks/07_integracion_cohete_y_riemann.ipynb` |
| 28/2-7/2 | Regla del trapezoide cohete | Encontrada | `notebooks/08_integracion_trapezoide_y_simpson.ipynb` |
| 28/2-7/2 | Regla del trapezoide compuesto cohete | Encontrada | `notebooks/08_integracion_trapezoide_y_simpson.ipynb` |
| 28/2-7/3 | Regla de simpson 1/3 cohete | Encontrada | `notebooks/08_integracion_trapezoide_y_simpson.ipynb` |
| 28/2-7/3 | Regla simpson 3/8 funcion propia | Encontrada | `notebooks/08_integracion_trapezoide_y_simpson.ipynb` |
| 14/3 | Metodo de la biseccion | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` |
| 14/3 | Metodo de Raphson | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` |
| 14/3 | Metodo de la secante | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` |
| 14/3-21/3 | Biseccion pelota flotante | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` |
| 14/3-21/3 | Newton-Raphson pelota flotante | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` |
| 14/3-21/3 | Secante pelota flotante | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` |
| 14/3-21/3 | Newton-Raphson ecuaciones simultaneas no lineales | Encontrada | `notebooks/09_raices_pelota_flotante.ipynb` |
| 28/3 | Ajuste por minimos cuadrados | No encontrada | No encontre notebook dedicado ni coincidencia clara. |
| 28/3 | Euler resolucion EDO | Encontrada | `notebooks/10_euler_ode_y_biseccion_decaimiento.ipynb` y `notebooks/14_metodo_disparo_y_odes_orden_superior.ipynb` |
| 28/3 | Regresion Exponencial | No encontrada | No hay coincidencia clara para regresion exponencial. |
| 28/3 | Regresion Lineal | Parcial/probable | Hay regresiones lineales en `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, pero no encontre una entrega de regresion lineal de marzo. |
| 28/3-4/4 | Biseccion para decaimiento | Encontrada | `notebooks/10_euler_ode_y_biseccion_decaimiento.ipynb` |
| 28/3-4/4 | ODE usando metodo de Euler | Encontrada | `notebooks/10_euler_ode_y_biseccion_decaimiento.ipynb` |
| 11/4 | Problema aleatorio | No encontrada | No aparece por nombre ni como seccion identificable. |
| 11/4 | RK2 | Encontrada | `notebooks/12_rk2_rk4_enfriamiento_y_burbuja.ipynb` |
| 11/4 | RK4 | Encontrada | `notebooks/11_runge_kutta_tareas.ipynb` y `notebooks/12_rk2_rk4_enfriamiento_y_burbuja.ipynb` |
| 11/4-18/4 | Tareas Runge Kutta | Encontrada | `notebooks/11_runge_kutta_tareas.ipynb` |
| 18/4-25/4 | Perfil de temperatura de la zona radiativa solar | Encontrada | `notebooks/13_diferencias_finitas_temperatura_solar.ipynb` |
| 18/4-2/4 | Pulsacion radial de una estrella | No encontrada | No encontre coincidencia clara. La fecha final `2/4/2026` parece inconsistente con el inicio `18/4/2026`. |
| 25/4 | Metodo de newton | No encontrada | Newton aparece en raices, pero no encontre una tarea de optimizacion con metodo de Newton. |
| 25/4 | Optimizacion canaleta | No encontrada | No encontre coincidencia clara. |
| 25/4 | Seccion aurea | Encontrada | `notebooks/15_seccion_aurea_temperatura_optima_emision_estelar.ipynb` |
| 25/4-2/5 | Temperatura optima de emision estelar | Encontrada | `notebooks/15_seccion_aurea_temperatura_optima_emision_estelar.ipynb` y tambien `17`. |
| 25/4-2/5 | Maximizar la luminosidad de un disco de acrecion | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, figura `17_disco_acrecion_gradiente.png`. |
| 25/4-2/5 | Problemas anteriores usando metodo de la gradiente | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb` implementa ascenso por gradiente. |
| 9/5 | Ecuacion difusion del calor | Encontrada | `notebooks/16_diferencias_finitas_enfriamiento_corteza_estrella_neutrones.ipynb` usa difusion termica FTCS. |
| 9/5 | Orbita leapfrog | No encontrada | No encontre coincidencia clara. |
| 9/5 | Orbita Verlet | No encontrada | No encontre coincidencia clara. |
| 9/5-16/5 | Enfriamiento termico de la corteza de una estrella de neutrones | Encontrada | `notebooks/16_diferencias_finitas_enfriamiento_corteza_estrella_neutrones.ipynb` |
| 16/5 | Cilindro de pared gruesa | No encontrada | No encontre coincidencia clara. |
| 16/5-23/5 | Potencial gravitacional en un disco protoplanetario | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, seccion Poisson disco. |
| 16/5-23/5 | Adveccion 1D | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, seccion MHD. |
| 16/5-23/5 | FEM para la ecuacion de Poisson gravitacional 1D | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, seccion FEM 1D. |
| 23/5-30/5 | Recreacion de ondas Alfven | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, figura `20_onda_alfven.png`. |
| 30/5 | Busqueda dentro de los datos de shows con SQL | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, funcion `task_sql`. |
| 30/5 | Busqueda dentro el SDSS | Parcial/probable | `17` menciona estructura de Astroquery/SDSS sin depender de red, pero no encontre resultados de consulta SDSS ejecutada. |
| 30/5-6/6 | Modificar los fits del notebook | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, salida `notebooks/17_tareas_pendientes_outputs/21_horsehead_modificado.fits`. |
| 30/5-6/6 | Busquedas en SDSS | Parcial/probable | Igual que SDSS diario: estructura presente, resultados no claros. |
| 30/5-6/6 | Los cuatro ejercicios de notebook | Parcial/probable | `17` cubre FITS, SQL y estructura SDSS; no pude verificar cuatro ejercicios separados. |
| 6/6 | Medidas estadisticas con numpy | Parcial/probable | `notebooks/05_todo_exercises01.ipynb` usa medias/estadisticas en datos de temperatura; no encontre tarea dedicada de ML. |
| 6/6 | Primer ejercicio de regresion lineal con ML | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, regresion Hubble. |
| 6/6-20/6 | Otros tres ejercicios de regresion lineal con ML | Parcial/probable | `17` contiene Hubble, SMBH y Bolshoi: tres regresiones en total; no pude confirmar cuatro ejercicios separados. |
| 13/6 | Clasificacion de tipos estelares | No encontrada | No encontre clasificador ni notebook ML 2 dedicado. |
| 13/6 | Clasificacion de galaxias, estrellas y QSO | Parcial/probable | `17` usa clases reales `GALAXY`, `QSO`, `STAR` para visualizar clustering, pero no implementa una clasificacion supervisada clara. |
| 13/6 | Relacion de los indices de color con el redshift | No encontrada | `17` calcula colores `u-g` y `g-z`, pero no encontre analisis con `redshift`. |
| 13/6 | Clasificacion estelar | No encontrada | No encontre clasificador estelar dedicado. |
| 27/6-25/7 | Clustering de galaxias con datos fotometricos reales | Encontrada | `notebooks/17_tareas_pendientes_mayo_julio.ipynb`, DBSCAN con datos fotometricos. |

## Tareas no encontradas claramente

- Ajuste por minimos cuadrados.
- Regresion exponencial.
- Problema aleatorio de EDO.
- Pulsacion radial de una estrella.
- Metodo de Newton para optimizacion.
- Optimizacion de canaleta.
- Orbita leapfrog.
- Orbita Verlet.
- Cilindro de pared gruesa.
- Clasificacion de tipos estelares.
- Relacion de indices de color con redshift.
- Clasificacion estelar.

## Parciales que conviene revisar

- `Exercise python / promedio de temperatura`: hay material de temperatura en `05` y preguntas de Exercises Python en `06`, pero el nombre exacto no aparece.
- `Regresion Lineal` de marzo: hay regresiones lineales en `17`, pero no encontre entrega separada de marzo.
- `SDSS`: `17` dice que deja estructura de Astroquery/SDSS sin red; no vi resultados concretos de busquedas SDSS.
- `Los cuatro ejercicios de notebook` de datos: hay FITS, SQL y SDSS parcial, pero no cuatro ejercicios claramente separados.
- `Otros tres ejercicios de regresion lineal con ML`: `17` tiene tres regresiones en total; si la lista esperaba cuatro ejercicios, falta confirmar uno.
- `Clasificacion galaxias/estrellas/QSO`: aparecen las clases reales en el clustering, pero no un modelo de clasificacion supervisada.

## Fechas sospechosas en la lista original

- `28/2/2026 -> 7/2/2026` para dos tareas de trapezoide: la fecha final es anterior a la inicial.
- `18/4/2026 -> 2/4/2026` para pulsacion radial: la fecha final es anterior a la inicial.
