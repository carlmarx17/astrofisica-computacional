# 📝 Tareas Resueltas — Astrofísica Computacional 2026-I

> Índice completo de las tareas del curso con resúmenes, métodos usados
> y enlaces directos a cada notebook.

---

## Resumen general

| #  | Tema | Métodos clave | Notebook |
|----|------|---------------|----------|
| 1  | Galaxia espiral + ejercicios numéricos | Modelado 2D, derivadas numéricas | [01](notebooks/01_galaxia_espiral_y_04_ejercicios.ipynb) |
| 2  | Altitud de Betelgeuse + repaso derivadas | Fórmula de altitud astronómica, diferencias centradas | [02](notebooks/02_clase_14_febrero_y_repaso.ipynb) |
| 3  | Diferencias finitas (entrega limpia) | Adelante, atrás, centradas; convergencia log-log | [03](notebooks/03_tarea_14_febrero_03_ejercicios.ipynb) |
| 4  | Cancelación catastrófica | Análisis de error de redondeo vs truncamiento | [04](notebooks/04_picos_verdes_en_grafica.ipynb) |
| 5  | TODOs de Exercises01 | Series de tiempo, clases, `.npz`/`.csv` | [05](notebooks/05_todo_exercises01.ipynb) |
| 6  | 5 preguntas de Exercises Python | Galaxia, Kepler, Altair, SunPy | [06](notebooks/06_respuestas_05_preguntas_exercises_python.ipynb) |
| 7  | Integración — cohete + Riemann | `scipy.integrate.quad`, sumas de Riemann | [07](notebooks/07_integracion_cohete_y_riemann.ipynb) |
| 8  | Trapezoide y Simpson | Trapezoide simple/compuesto, Simpson 1/3 y 3/8 | [08](notebooks/08_integracion_trapezoide_y_simpson.ipynb) |
| 9  | Búsqueda de raíces — pelota flotante | Bisección, Newton-Raphson, Secante | [09](notebooks/09_raices_pelota_flotante.ipynb) |
| 10 | Euler + bisección para decaimiento | Euler explícito, decaimiento C-14 | [10](notebooks/10_euler_ode_y_biseccion_decaimiento.ipynb) |
| 11 | Runge-Kutta — 4 problemas | RK4, Secante, Bisección, Newton-Raphson | [11](notebooks/11_runge_kutta_tareas.ipynb) |
| 12 | RK2/RK4 — enfriamiento + burbuja | Ralston, Runge, Kutta; bisección | [12](notebooks/12_rk2_rk4_enfriamiento_y_burbuja.ipynb) |
| 13 | Perfil Temp. Solar | Diferencias finitas, convergencia | [13](notebooks/13_diferencias_finitas_temperatura_solar.ipynb) |
| 14 | ODEs Orden Superior | Método de disparo, reducción de orden | [14](notebooks/14_metodo_disparo_y_odes_orden_superior.ipynb) |
| 15 | Temperatura óptima de emisión estelar | Sección áurea, optimización unidimensional | [15](notebooks/15_seccion_aurea_temperatura_optima_emision_estelar.ipynb) |
| 16 | Enfriamiento de corteza de estrella de neutrones | FTCS, difusión térmica, pérdidas por neutrinos | [16](notebooks/16_diferencias_finitas_enfriamiento_corteza_estrella_neutrones.ipynb) |
| 17 | Pendientes mayo-julio | Optimización, EDP, FEM, MHD, FITS, SQL, ML, clustering | [17](notebooks/17_tareas_pendientes_mayo_julio.ipynb) |

---

## Detalle por tarea

### 1 · Galaxia espiral y 04 ejercicios

**Tema:** Modelado visual · Derivadas numéricas  
**Notebook:** [`01_galaxia_espiral_y_04_ejercicios.ipynb`](notebooks/01_galaxia_espiral_y_04_ejercicios.ipynb)

Se construyó un modelo visual de galaxia espiral con bulbo central, disco, barra y
modulación de brazos espirales. Aparte se resolvieron los tres ejercicios numéricos
de `04. Ejercicios.ipynb`.

| Resultado | Valor |
|-----------|-------|
| $a_0$ | 10.4 km s⁻² |
| $a_4$ | 13.2 km s⁻² |
| $a_2$ | 11.7 km s⁻² |
| $dP/dr$ | −12.9 |
| $dz/dD$ | 0.001125 Mpc⁻¹ |
| $d^2z/dD^2$ | 0.000125 Mpc⁻² |

---

### 2 · Altitud de Betelgeuse y repaso de derivadas

**Tema:** Astronomía de posición · Diferencias centradas  
**Notebook:** [`02_clase_14_febrero_y_repaso.ipynb`](notebooks/02_clase_14_febrero_y_repaso.ipynb)

**Parte A — Gráfica de clase.** Altitud de Betelgeuse usando:

$$\sin h = \sin\delta\,\sin\phi + \cos\delta\,\cos\phi\,\cos H$$

Se barrió el ángulo horario $H \in [0, 2\pi]$ para dos ciudades.

**Parte B — Segunda derivada.** Diferencias centradas $\mathcal{O}(h^2)$:

$$f''(x) \approx \frac{f(x+h) - 2f(x) + f(x-h)}{h^2}$$

| Resultado | Bogotá | Medellín |
|-----------|--------|----------|
| $h_{\max}$ | 87.25° | 88.85° |
| Tiempo sobre 30° | ~8.04 h | ~8.06 h |

El error de la segunda derivada escala como $\mathcal{O}(h^2)$, confirmado numéricamente.

---

### 3 · Diferencias finitas — entrega limpia

**Tema:** Derivadas numéricas · Convergencia  
**Notebook:** [`03_tarea_14_febrero_03_ejercicios.ipynb`](notebooks/03_tarea_14_febrero_03_ejercicios.ipynb)

Versión final (entrega) de los tres ejercicios: diferencias hacia adelante, hacia atrás
y centradas para primera y segunda derivada. Se graficó error vs $h$ en escala log-log
para verificar el orden de convergencia.

| Resultado | Valor |
|-----------|-------|
| $a_0$ | 10.4 km s⁻² |
| $a_4$ | 13.2 km s⁻² |
| $a_2$ | 11.7 km s⁻² |
| $dP/dr$ | −12.9 |

---

### 4 · Picos verdes — cancelación catastrófica

**Tema:** Aritmética de punto flotante  
**Notebook:** [`04_picos_verdes_en_grafica.ipynb`](notebooks/04_picos_verdes_en_grafica.ipynb)

Explicación de por qué la curva de error (verde) de diferencias finitas se dispara
para $h$ muy pequeños: la resta $f(x+h) - f(x)$ pierde dígitos significativos
(**cancelación catastrófica**) y el error de redondeo domina sobre el de truncamiento.

---

### 5 · TODOs de Exercises01

**Tema:** Análisis de datos · Clases en Python  
**Notebook:** [`05_todo_exercises01.ipynb`](notebooks/05_todo_exercises01.ipynb)

Se completaron todos los `# TODO` del notebook `Exercises01.ipynb`:

- Filtrado de outliers en datos de temperatura
- Cálculo de anomalías con media móvil
- Clases `TimeSeriesAnalyzer` y `WeatherAnalyzer`
- Exportación a `.npz` y `.csv`

| Dato | Valor |
|------|-------|
| Outliers removidos | 2 (fracción: 0.00548) |
| Reducción de $\sigma$ con suavizado | ~6.73 % |

---

### 6 · Respuestas a las 5 preguntas de Exercises Python

**Tema:** Galaxia, Kepler, Altair, SunPy  
**Notebook:** [`06_respuestas_05_preguntas_exercises_python.ipynb`](notebooks/06_respuestas_05_preguntas_exercises_python.ipynb)

Se respondieron las 5 preguntas abiertas: galaxia espiral (`contourf` logarítmico),
altitud de Betelgeuse, ajuste de órbita de Kepler con `least_squares`,
visualización interactiva en Altair, y análisis de oscilaciones coronales con SunPy.

---

### 7 · Integración — cohete analítico y sumas de Riemann

**Tema:** Integración numérica  
**Notebook:** [`07_integracion_cohete_y_riemann.ipynb`](notebooks/07_integracion_cohete_y_riemann.ipynb)  
**Origen:** `01. Fundamental Algorithms / 03. Integration / Integracion.pdf`

- **Cohete:** distancia entre $t=8$ y $t=30$ s con `scipy.integrate.quad`.
- **Riemann:** sumas a la derecha para $y = x^3$ en $[0, b]$, convergencia $\mathcal{O}(1/n)$.

| Resultado | Valor |
|-----------|-------|
| Distancia cohete (exacta) | ≈ 11 061.24 m |
| Riemann $b=4$, $n=1000$ | error < 0.001 |

---

### 8 · Trapezoide y Simpson para el cohete

**Tema:** Integración numérica avanzada  
**Notebook:** [`08_integracion_trapezoide_y_simpson.ipynb`](notebooks/08_integracion_trapezoide_y_simpson.ipynb)  
**Origen:** `01. Fundamental Algorithms / 03. Integration / Integracion.pdf`

Se aplicaron cuatro métodos a $\int_8^{30} v(t)\,dt$ y a una función propia:

| Método | Aproximación | Error rel. |
|--------|-------------|------------|
| Trapezoide simple | alta sobreestimación | ~15 % |
| Trapezoide compuesto ($n=2$) | mucho mejor | ~4 % |
| Simpson 1/3 simple | muy preciso | < 0.01 % |
| Simpson 3/8 (f propia) | exacto a 6 cifras | < 0.001 % |

---

### 9 · Búsqueda de raíces — pelota flotante

**Tema:** Métodos de búsqueda de raíces  
**Notebook:** [`09_raices_pelota_flotante.ipynb`](notebooks/09_raices_pelota_flotante.ipynb)  
**Origen:** `01. Fundamental Algorithms / 04. RootSearching / Roots_searching.pdf`

Ecuación: $f(x) = x^3 - 0.165\,x^2 + 3.993 \times 10^{-4} = 0$

| Método | Detalle | Resultado |
|--------|---------|-----------|
| Bisección | $[0, 0.05]$ | $x \approx 0.02938$ m |
| Newton-Raphson | $x_0 = 0.05$ | converge en 4 it. |
| Secante | $x_0=0,\; x_1=0.05$ | converge en ~5 it. |
| NR (sistema 2×2) | $xy=2,\; x^2+y=5$ | $(x,y) \approx (1.0,\, 2.0)$ |

---

### 10 · Euler y bisección para decaimiento radiactivo

**Tema:** EDOs · Decaimiento C-14  
**Notebook:** [`10_euler_ode_y_biseccion_decaimiento.ipynb`](notebooks/10_euler_ode_y_biseccion_decaimiento.ipynb)  
**Origen:** `02. DifferentialEquations / EulerMethod / odes_and_euler.pdf`

1. **Euler:** $3\,dy/dx + 2y = e^{-x}$, $y(0)=5$. Solución exacta $y = 6e^{-2x/3} - e^{-x}$.
   Convergencia $\mathcal{O}(h)$ verificada.
2. **Bisección C-14:** hallar $t^*$ tal que $N(t^*)/N_0 = 10\%$, con $T_{1/2} = 5730$ años.

| Resultado | Valor |
|-----------|-------|
| Error máximo Euler ($h=0.1$) | ≈ 0.036 |
| $t^*$ decaimiento C-14 | ≈ 19 035 años |

---

### 11 · Runge-Kutta — 4 problemas aplicados

**Tema:** RK4 + búsqueda de raíces  
**Notebook:** [`11_runge_kutta_tareas.ipynb`](notebooks/11_runge_kutta_tareas.ipynb)  
**Origen:** `02. DifferentialEquations / Runge-Kutta / runge_kutta.pdf`

| # | Problema | ODE | Método raíces | Resultado |
|---|----------|-----|---------------|-----------|
| 1 | Circuito RC — media carga | $dq/dt = (V - q/C)/R$ | Secante | $t_m = RC\ln 2 \approx 0.693$ s |
| 2 | Transferencia radiativa | $dI/d\tau = -I + S$ | Bisección | $\tau_m \approx 0.847$ |
| 3 | Fulguración solar | $dE/dt = -\alpha E^n$ | Newton-Raphson | $t_m \approx 0.828$ s |
| 4 | Burbuja ascendente — 99 % $v_t$ | $dv/dt = A - Bv$ | Bisección | $t_m \approx 0.00123$ s |

Todos validados contra la solución analítica.

---

### 12 · RK2 y RK4 — enfriamiento radiativo y burbuja ascendente

**Tema:** Variantes de Runge-Kutta · Bisección  
**Notebook:** [`12_rk2_rk4_enfriamiento_y_burbuja.ipynb`](notebooks/12_rk2_rk4_enfriamiento_y_burbuja.ipynb)

#### Parte A — Enfriamiento radiativo

EDO: $\;d\theta/dt = -A\,(\theta^4 - B)$, con $\theta_0 = 1200$ K y $h = 60$ s.

Se compararon tres integradores:

| Método | $\theta(480\;\text{s})$ |
|--------|------------------------|
| RK2 (Ralston) | 652.2548 K |
| RK4 (Runge clásico) | 647.5393 K |
| RK4 (Kutta) | 647.4437 K |

#### Parte B — Burbuja ascendente (detalle)

EDO: $\;dv/dt = A - Bv$ (Stokes). Se integró con RK4 y se usó bisección para
localizar el instante exacto en que $v = 0.99\,v_t$.

| Resultado | Valor |
|-----------|-------|
| $v_t$ | 2.17738 m/s |
| $t_m$ (bisección, 19 iteraciones) | 0.001228 s |

---

### 13 · Perfil de Temperatura Solar (Diferencias Finitas)

**Tema:** Diferencias finitas para BVP  
**Notebook:** [`13_diferencias_finitas_temperatura_solar.ipynb`](notebooks/13_diferencias_finitas_temperatura_solar.ipynb)

Resolución de la ecuación de Laplace radial para la temperatura interior del Sol usando diferencias finitas. Se analiza la convergencia para determinar el número de nodos $N$ requerido para un error relativo $<1\%$ y $<0.01\%$.

---

### 14 · Método de Disparo y ODEs de Orden Superior

**Tema:** Método de disparo · Reducción de orden  
**Notebook:** [`14_metodo_disparo_y_odes_orden_superior.ipynb`](notebooks/14_metodo_disparo_y_odes_orden_superior.ipynb)

---

### 15 · Temperatura óptima de emisión estelar

**Tema:** Optimización numérica · Sección áurea  
**Notebook:** [`15_seccion_aurea_temperatura_optima_emision_estelar.ipynb`](notebooks/15_seccion_aurea_temperatura_optima_emision_estelar.ipynb)

Se maximiza la función

$$
P(T) = \sigma T^4 e^{-T/T_0}\left(1 - e^{-h\nu/(k_B T)}\right)
$$

en el intervalo $T \in [3000, 50000]$ K usando el método de la **sección áurea**
con tolerancia $\varepsilon = 50$ K. Además, se grafica la curva $P(T)$ y se
marca la temperatura óptima hallada.

| Resultado | Valor |
|-----------|-------|
| Iteraciones | 15 |
| $T_\mathrm{opt}$ | $\approx 3.08 \times 10^4$ K |
| $P(T_\mathrm{opt})$ | $\approx 3.51 \times 10^8$ W m$^{-2}$ |

Dos ejemplos aplicados de clase:
1. Resolución de un sistema derivado de una ODE de orden superior usando una aproximación inicial para $z(13)$ y avance con método de Euler.
2. Reducción de $3y'' + 2y' + 5y = e^{-x}$ a un sistema de primer orden, resolviéndolo en Python y comparando los resultados de los métodos de Euler y Runge-Kutta de 4to orden (RK4).

---

### 16 · Enfriamiento térmico de corteza de estrella de neutrones

**Tema:** EDP parabólicas · FTCS  
**Notebook:** [`16_diferencias_finitas_enfriamiento_corteza_estrella_neutrones.ipynb`](notebooks/16_diferencias_finitas_enfriamiento_corteza_estrella_neutrones.ipynb)  
**Script:** [`16_diferencias_finitas_enfriamiento_corteza_estrella_neutrones.py`](notebooks/16_diferencias_finitas_enfriamiento_corteza_estrella_neutrones.py)

Se resuelve la difusión térmica 1D en una corteza de estrella de neutrones:

$$
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial z^2} - \epsilon_\nu T^5
$$

con fronteras de Dirichlet, esquema explícito FTCS y término de enfriamiento por neutrinos.

---

### 17 · Pendientes mayo-julio 2026

**Tema:** Optimización · EDP · FEM · MHD · FITS · SQL · ML  
**Notebook:** [`17_tareas_pendientes_mayo_julio.ipynb`](notebooks/17_tareas_pendientes_mayo_julio.ipynb)  
**Script:** [`17_tareas_pendientes_mayo_julio.py`](notebooks/17_tareas_pendientes_mayo_julio.py)  
**Resumen:** [`17_tareas_pendientes_mayo_julio_resumen.md`](notebooks/17_tareas_pendientes_mayo_julio_resumen.md)

Entrega consolidada para las tareas pendientes detectadas en la lista del profesor:

- Luminosidad de disco de acreción y método del gradiente.
- Potencial gravitacional en disco protoplanetario con Gauss-Seidel y SOR.
- FEM 1D para Poisson gravitacional.
- Advección 1D de campo magnético y onda de Alfvén.
- Modificación de FITS local, ejercicios SQL, regresiones de ML y clustering fotométrico con DBSCAN.

| Resultado | Valor |
|-----------|-------|
| $T_\mathrm{opt}$ | $\approx 30794$ K |
| SOR disco protoplanetario | 79 iteraciones |
| DBSCAN | 3 clusters, 5302 puntos de ruido |

---

*Última actualización: julio 2026*
