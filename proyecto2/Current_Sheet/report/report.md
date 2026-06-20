# Proyecto 2: Magnetohidrodinámica - Simulación y análisis
## Harris Current Sheet con Hall MHD

---

## Resumen y guía de lectura

Este reporte documenta una simulación de reconexión magnética en una lámina de corriente de Harris usando Hall MHD con PLUTO. La parte de PLUTO resuelve la evolución física completa; la parte de Python reproduce la condición inicial, calcula diagnósticos y genera las figuras usadas para interpretar la corrida.

El proyecto queda organizado así:

| Ruta | Contenido |
|------|-----------|
| `Current_Sheet/` | Caso principal: configuración PLUTO, análisis, figuras y reporte. |
| `Current_Sheet/pluto_sim/` | Copia de la configuración usada en la corrida y, si están presentes localmente, salidas VTK/DBL regenerables. |
| `Current_Sheet/analysis/` | Scripts de postproceso para leer PLUTO, calcular diagnósticos y producir figuras finales. |
| `Current_Sheet/python_reproduction/` | Reproducción Python de la condición inicial y herramientas de chequeo. |
| `Whistler_Waves/` | Configuraciones auxiliares para ondas whistler. |
| `PLUTO/` | Código PLUTO usado como dependencia externa/local. |

Los resultados centrales son: flujo reconectado final $4.5525$ bajo la métrica $\int_0^{L_x/2}|B_y(x,0)|dx$, máximo temporal $\max |J_z|=3.0819$ cerca de $t\approx25$, y error final $||\nabla\cdot\mathbf{B}||_2=4.01\times10^{-4}$. La comparación Python vs PLUTO en $t=0$ valida que la condición inicial fue reproducida con errores relativos L2 de orden $10^{-8}$.

### Correspondencia con el enunciado

| Requisito del proyecto | Dónde se responde | Estado |
|------------------------|-------------------|--------|
| Elegir un problema MHD al menos 2D desde `PLUTO/Test_Problems/MHD` | Se usa `MHD/Hall_MHD/Current_Sheet`, una lámina de Harris 2D. | Cumplido |
| Marco teórico, ecuaciones, condiciones de frontera y contexto físico | Secciones 1.1-1.5. | Cumplido |
| Reproducir una simulación PLUTO completa | Secciones 2.1-2.6 y archivos `Current_Sheet/init.c`, `definitions_01.h`, `pluto_01.ini`. | Cumplido |
| Graficar variables físicas con pyPLUTO | Secciones 2.3, 2.5 y script `analysis/plot_results.py`. | Cumplido |
| Implementación Python independiente | `python_reproduction/hall_mhd_harris.py` reproduce la condición inicial y el postproceso; no implementa un solver Hall-MHD completo. | Parcial |
| Comparación PLUTO vs Python y errores | Secciones 4.1-4.4 y `analysis/analysis.py`. | Cumplido para la condición inicial y diagnósticos |
| Discusión de limitaciones, fuentes de error y mejoras | Sección 5. | Cumplido |

La principal diferencia frente al enunciado ideal es que la parte Python no reemplaza a PLUTO como solver evolutivo Hall-MHD. En este trabajo Python se usa como reproducción independiente del setup y como herramienta de análisis cuantitativo. Un solver Hall-MHD completo requiere un esquema conservativo 2D con Riemann solver, control estricto de $\nabla\cdot\mathbf{B}$ y tratamiento estable del término Hall; por eso se declara explícitamente como mejora futura.

---

## 1. Marco Teórico

### 1.1 Ecuaciones de la MHD Hall

La magnetohidrodinámica (MHD) describe un plasma como un fluido conductor gobernado por las ecuaciones de conservación de masa y momento acopladas a las ecuaciones de Maxwell. En el formalismo **MHD Hall**, se incluye el término de corriente de Hall en la ley de inducción, que separa el movimiento del fluido electrónico del iónico.

Las ecuaciones MHD Hall isotérmicas en 2D son:

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0
$$

$$
\frac{\partial (\rho \mathbf{v})}{\partial t} + \nabla \cdot \left( \rho \mathbf{v} \mathbf{v} + P_t \mathbb{I} - \mathbf{B} \mathbf{B} \right) = 0
$$

$$
\frac{\partial \mathbf{B}}{\partial t} - \nabla \times \left( \mathbf{v} \times \mathbf{B} \right) + \nabla \times \left( \frac{\mathbf{J} \times \mathbf{B}}{n_e e} \right) = 0
$$

donde $P_t = P + B^2/2$ es la presión total, $P = c_s^2 \rho$ (EOS isotérmica), $\mathbf{J} = \nabla \times \mathbf{B}$ es la densidad de corriente, y el término $\mathbf{J} \times \mathbf{B} / (n_e e)$ es el término Hall.

La ecuación de inducción puede reescribirse como:

$$
\frac{\partial \mathbf{B}}{\partial t} = \nabla \times \left( \mathbf{v}_e \times \mathbf{B} \right), \quad \mathbf{v}_e = \mathbf{v} - \frac{\mathbf{J}}{n_e e}
$$

donde $\mathbf{v}_e$ es la velocidad del fluido electrónico. El término Hall congela las líneas de campo en el fluido electrónico en lugar del iónico, permitiendo fenómenos como la **reconexión magnética rápida** y la propagación de **ondas whistler**.

### 1.2 La Lámina de Harris

La configuración de Harris (Harris, 1962) es un equilibrio MHD unidimensional que consiste en una lámina de corriente donde el campo magnético invierte su dirección:

$$
\mathbf{B}(y) = B_0 \tanh\left(\frac{y}{l}\right) \hat{\mathbf{x}}, \quad \rho(y) = \rho_0 + \frac{\rho_1}{\cosh^2(y/l)}
$$

Para desencadenar la reconexión, se superpone una perturbación magnética de la forma:

$$
\delta B_x = -\Psi_0 k_y \sin(k_y y) \cos(2k_x x), \quad \delta B_y = 2\Psi_0 k_x \sin(2k_x x) \cos(k_y y)
$$


### 1.3 El Efecto Hall en la Reconexión

El término Hall juega un papel crucial en la reconexión magnética:
- Permite la separación de escalas iónicas y electrónicas
- Genera ondas whistler que transportan información rápidamente
- Acelera la reconexión comparado con MHD ideal/resistivo
- Crea estructuras coherentes en la región de difusión

### 1.4 Condiciones de frontera y normalización

La corrida se hace en geometría cartesiana 2D con coordenadas $(x,y)$ y una tercera dirección inactiva. El dominio físico es

$$
x\in[-12.8,12.8], \qquad y\in[-6.4,6.4],
$$

discretizado con una malla uniforme $256\times128$. Las fronteras usadas en `pluto_01.ini` son periódicas en $x$ y reflectivas en $y$. La periodicidad en $x$ permite que la perturbación inicial sea compatible con el dominio horizontal; las fronteras reflectivas en $y$ mantienen confinada la lámina dentro del dominio vertical. La dirección $z$ se deja con una sola celda y condiciones de salida porque el problema es estrictamente 2D.

Las variables están en unidades normalizadas de PLUTO. Se usa una ecuación de estado isotérmica con $c_s^2=0.5$, densidad de fondo $\rho_{\rm bg}=0.2$, campo característico $B_0=1$, ancho de lámina $l=0.5$ y amplitud de perturbación $\Psi_0=0.02$.

### 1.5 Estado del arte resumido

La lámina de Harris es una configuración estándar para estudiar reconexión porque contiene una inversión de campo magnético sostenida por una capa de corriente. Harris (1962) introdujo este equilibrio como modelo idealizado de una hoja de plasma. En reconexión magnética moderna, el problema se usa para estudiar cómo cambia la topología del campo y cómo se convierte energía magnética en energía cinética y térmica.

El desafío GEM de Birn et al. (2001) convirtió la reconexión en una prueba comparativa entre distintos modelos numéricos: MHD resistiva, Hall MHD, híbridos y cinéticos. Una conclusión central de esa línea de trabajo es que el término Hall permite reconexión más rápida que la MHD resistiva simple, porque desacopla el movimiento electrónico del iónico cerca de la región de difusión. Huba (2003) resume la física Hall y muestra la conexión con ondas whistler, que transportan información magnética a escalas pequeñas. PLUTO incorpora estos ingredientes en un marco conservativo de dinámica de fluidos astrofísicos (Mignone et al. 2012), lo que permite reproducir pruebas 2D como la lámina de Harris con configuraciones controladas.

En aplicaciones astrofísicas, Hall MHD aparece en plasmas parcialmente ionizados, discos protoplanetarios, magnetosferas y evolución magnética de objetos compactos. Trabajos como Lesur et al. (2014) y Viganò et al. (2012) muestran que el término Hall puede modificar la estabilidad, el transporte angular y la evolución del campo magnético. Por eso, aunque este proyecto usa una prueba idealizada, el mecanismo físico estudiado es relevante para problemas más generales de reconexión y transporte magnético.

---

## 2. Simulación PLUTO

### 2.1 Configuración

| Parámetro | Valor |
|-----------|-------|
| Código | PLUTO v2026 |
| Física | MHD Hall isotérmico |
| Dimensiones | 2D Cartesiano |
| Dominio | $[-12.8, 12.8] \times [-6.4, 6.4]$ |
| Grilla | $256 \times 128$ |
| EOS | Isotérmica ($c_s^2 = 0.5$) |
| Solver | HLL |
| Time stepping | RK2 |
| Div(B) control | DIV_CLEANING |
| Hall MHD | EXPLICIT |
| CFL | 0.25 |
| $l$ (ancho lámina) | 0.5 |
| $\Psi_0$ (perturbación) | 0.02 |
| $t_{stop}$ | 60.0 |

### 2.2 Resultados

La simulación evoluciona la lámina de Harris desde el equilibrio perturbado hasta la reconexión completa:

- **t = 0**: Estado inicial con la lámina de corriente y la perturbación senoidal
- **t ≈ 15-20**: Comienza la reconexión; las líneas de campo se deforman
- **t ≈ 30-40**: Reconexión activa; se forman islas magnéticas
- **t ≈ 50-60**: Estado cuasi-estacionario con flujo reconectado significativo

El diagnostico principal se calculo como el flujo reconectado sobre el eje medio, $\int_0^{L_x/2}|B_y(x,0)|dx$. La curva en `analysis/figs/diagnostics_timeseries.png` muestra crecimiento claro desde $t\approx15$ y alcanza $4.55$ en $t=60$. La corriente maxima $|J_z|$ aumenta durante el onset de reconexion, con maximo cercano a $3.08$ en $t\approx25$, y despues decrece cuando la configuracion entra en una fase mas relajada.

### 2.3 Visualización con pyPLUTO

Se generaron 13 snapshots cubriendo $t=0$ a $t=60$ con todas las variables fisicas ($\rho$, $P$, $v_x$, $v_y$, $B_x$, $B_y$, $B_z$, $|\mathbf{B}|$). El panel `analysis/figs/pluto_final_all_variables.png` resume el estado final e incluye tambien $J_z=\partial_x B_y-\partial_y B_x$. Como la salida VTK se guardo cada 5 unidades de tiempo, el snapshot mas cercano a $t=57$ es $t=55$; la figura `analysis/figs/jz_fieldlines_t55.png` muestra $J_z$ con lineas de campo magnetico y anotaciones de los puntos O/X.

![Diagnosticos temporales: flujo reconectado, corriente maxima y divergencia de B](../analysis/figs/diagnostics_timeseries.png)

**Figura 1.** Diagnosticos temporales de la corrida Hall MHD. El flujo reconectado crece despues de $t\approx20$, la corriente maxima $|J_z|$ alcanza su mayor valor cerca de $t\approx25$, y el error solenoidal $||\nabla\cdot\mathbf{B}||_2$ permanece acotado durante toda la evolucion.

![Estado final con todas las variables fisicas](../analysis/figs/pluto_final_all_variables.png)

**Figura 2.** Estado final en $t=60$ para densidad, presion, velocidades, campo magnetico, modulo de campo y corriente $J_z$.

![Corriente y lineas de campo](../analysis/figs/jz_fieldlines_t55.png)

**Figura 3.** Densidad de corriente $J_z$ con lineas de campo magnetico. El marcador O indica el centro aproximado de la isla magnetica y los marcadores X indican regiones laterales donde cambia la conectividad del campo. El snapshot mas cercano al tiempo de referencia $t=57$ es $t=55$ por la cadencia de salida.

![Firma Hall en t=60](../analysis/figs/hall_signature_t60.png)

**Figura 4.** Diagnosticos especificos del efecto Hall en $t=60$: campo fuera del plano $B_z$, componente fuera del plano del termino Hall $(\mathbf{J}\times\mathbf{B})_z/\rho$ y magnitud del desacoplamiento ion-electron $|\mathbf{v}_e-\mathbf{v}|=|\mathbf{J}|/\rho$, usando $n_e e\simeq\rho$ en las unidades normalizadas de esta corrida.

### 2.4 Interpretacion de los diagnosticos

Las salidas muestran de forma clara el proceso de reconexion magnetica en una lamina de Harris con Hall MHD. La simulacion empieza cerca de un equilibrio perturbado, luego la lamina se deforma, aparece reconexion y finalmente se alcanza una fase no lineal con islas magneticas.

**Flujo reconectado.** La primera grafica de la Figura 1 muestra que el flujo reconectado crece muy poco al inicio y luego aumenta rapidamente entre $t\approx25$ y $t\approx45$. Esto indica una fase de crecimiento fuerte de la reconexion. Despues de $t\approx45$, el flujo se estabiliza alrededor de $\psi_{\rm rec}\approx4.4-4.8$. Fisicamente, la evolucion puede resumirse como:

$$
\text{inicio lento} \rightarrow \text{reconexion rapida} \rightarrow \text{saturacion no lineal}.
$$

Este comportamiento es razonable para una lamina de corriente que se rompe y forma islas magneticas.

**Intensidad de la lamina de corriente.** La segunda grafica muestra $\max |J_z|$. La corriente aumenta hasta un maximo cercano a $t\approx25$, con $\max |J_z|\approx3.1$, y luego cae de forma marcada despues de $t\approx35-45$. Esto tiene sentido fisico: al inicio la lamina se comprime y la corriente se intensifica; despues, cuando la reconexion ya esta desarrollada, la lamina original se rompe, se ensancha o se reorganiza en estructuras tipo islas/plasmoides. Por lo tanto, la caida de $\max |J_z|$ no indica un fallo numerico por si misma; puede indicar que la lamina inicial ya fue destruida por la reconexion.

**Error de divergencia.** La tercera grafica muestra $||\nabla\cdot\mathbf{B}||_2$. El error parte de valores muy pequenos, sube hasta el orden de $10^{-3}$ durante la fase no lineal y luego baja hacia valores de orden $4\times10^{-4}-6\times10^{-4}$. En MHD idealmente se exige $\nabla\cdot\mathbf{B}=0$, por lo que esta cantidad debe reportarse. En esta corrida el error no parece catastrofico: crece cuando la dinamica se vuelve no lineal, pero se mantiene relativamente controlado. Como se usa limpieza de divergencia, este comportamiento es esperable; con un esquema estrictamente solenoidal tipo constrained transport se esperarian errores aun mas bajos.

**Mapa de $J_z$ y lineas de campo en $t=55$.** La Figura 3 es la evidencia visual mas importante. Se observa una isla magnetica central alrededor de $x\approx0$, $y\approx0$: las lineas de campo se cierran alrededor de esa region, firma de una estructura tipo O-point o plasmoide. Tambien se ven regiones tipo X-point aproximadamente a los lados de la isla central, cerca de $x\approx-5$ y $x\approx5$, donde las lineas de campo cambian de conectividad. Esto es precisamente lo esperado en reconexion magnetica. Las acumulaciones fuertes de $J_z$ en la isla central y hacia los bordes laterales pueden estar asociadas tanto a la periodicidad del dominio como a la formacion de islas adicionales en los extremos.

**Firma Hall.** La Figura 4 muestra una huella mas directa del efecto Hall. Aunque la condicion inicial tiene $B_z=0$ y $v_z=0$, la evolucion genera campo y velocidad fuera del plano: en $t=60$ se obtiene $\max |B_z|=0.223$ y $\max |v_z|=0.445$. La componente $(\mathbf{J}\times\mathbf{B})_z/\rho$ alcanza valores de orden $0.208$, lo que indica que el termino Hall no solo esta activado numericamente, sino que contribuye de forma localizada cerca de la region de reconexion.

La secuencia fisica global es:

$$
\text{lamina de corriente} \rightarrow \text{intensificacion de } J_z \rightarrow \text{ruptura de la lamina} \rightarrow \text{X-points e islas magneticas} \rightarrow \text{saturacion del flujo reconectado}.
$$

Un detalle importante es que la cantidad llamada flujo reconectado en esta entrega se calcula como

$$
\int_0^{L_x/2}|B_y(x,0)|\,dx,
$$

lo cual funciona como indicador global de reconexion, pero no es exactamente el flujo reconectado clasico. Para obtener una medida mas local se reconstruyo tambien el potencial vectorial $A_z$ a partir de

$$
B_x=\frac{\partial A_z}{\partial y}, \qquad B_y=-\frac{\partial A_z}{\partial x},
$$

y se estimo

$$
\psi_{A_z}=|A_z(O)-A_z(X)|.
$$

En $t=60$, usando el O-point central y el X-point lateral mas contrastado, se obtiene $\psi_{A_z}\approx0.915$. Esta cantidad no reemplaza por completo un algoritmo automatico robusto de deteccion de puntos criticos, pero es una medida fisicamente mas cercana al flujo reconectado clasico que la integral global de $|B_y|$.

### 2.5 Evolucion fisica por tiempos

La evolucion visual es coherente con una lamina de corriente de Harris en Hall MHD: inicia en equilibrio, la perturbacion deforma la lamina, aparece una region tipo X-point, crece la reconexion y finalmente se desarrolla una fase no lineal con estructuras magneticas tipo islas/plasmoides.

| Tiempo | Estado fisico | Evidencia principal |
|--------|---------------|--------------------|
| $t=0$ | Equilibrio perturbado | Lamina de Harris, $v\approx0$, inversion de $B_x$ y perturbacion inicial en $B_y$. |
| $t=20$ | Inicio de reconexion | Deformacion de la lamina, aparicion de flujos y crecimiento de la componente reconectada. |
| $t=30$ | Reconexion activa | Intensificacion de $J_z$, reorganizacion de la capa de corriente y region tipo X. |
| $t=50$ | Regimen no lineal | Isla magnetica, redistribucion de densidad y estructura compleja en $J_z$. |
| $t=60$ | Saturacion aproximada | Flujo reconectado casi estabilizado y presencia clara de componentes Hall fuera del plano. |

![Snapshot t=0](../plots/hall_cs_0000_t0.0.png)

**$t=0$.** La condicion inicial esta bien representada. La densidad $\rho$ y la presion $P=c_s^2\rho$ se concentran alrededor de $y=0$, donde esta la lamina de corriente. El campo $B_x$ cambia de signo al cruzar $y=0$, como corresponde a una configuracion de Harris. La corriente $J_z$ aparece concentrada en la lamina porque $B_x$ varia bruscamente en la direccion $y$. Las velocidades $v_x$ y $v_y$ son practicamente nulas, y $B_y$ contiene solo la perturbacion inicial que sirve como semilla de reconexion.

![Snapshot t=20](../plots/hall_cs_0004_t20.0.png)

**$t=20$.** La reconexion empieza a ser visible. La lamina de densidad y presion se curva, aparecen velocidades en ambas direcciones y $B_y$ deja de ser una perturbacion puramente inicial para mostrar estructura de campo reconectado. La corriente $J_z$ todavia conserva una banda elongada sobre $y=0$, pero ya se modifica cerca del centro del dominio, indicando la formacion de una region tipo X-point.

![Snapshot t=30](../plots/hall_cs_0006_t30.0.png)

**$t=30$.** La reconexion esta mas desarrollada. La lamina se estrecha en la zona central, la densidad y presion se redistribuyen hacia regiones de salida, y $|B|$ muestra zonas de campo reducido alrededor de la region central. La corriente $J_z$ deja de ser una banda uniforme y se concentra en estructuras mas localizadas, senal de reorganizacion de la capa de corriente.

![Snapshot t=50](../plots/hall_cs_0010_t50.0.png)

**$t=50$.** El sistema entra en una fase no lineal. La densidad ya no aparece como una lamina continua sino como acumulaciones localizadas; $B_y$ y $|B|$ muestran estructuras cerradas o tipo isla magnetica; y $J_z$ adquiere una geometria compleja, compatible con regiones tipo X y O. Las velocidades tambien son mas estructuradas, lo que indica que el plasma ya no responde como una perturbacion lineal simple.

### 2.6 Costo computacional y salida de datos

La corrida PLUTO fue ejecutada en un solo procesador sobre una malla uniforme $256\times128$. El log de ejecucion registra los siguientes datos:

| Cantidad | Valor |
|----------|------:|
| Tiempo fisico simulado | $t=60$ |
| Pasos hidrodinamicos/MHD | 194828 |
| Snapshots VTK | 13 |
| Intervalo de salida VTK | $\Delta t=5$ |
| Memoria asignada | 11.22 MB |
| Tiempo de ejecucion medido por PLUTO | 1 h 22 min 39 s |
| Tiempo promedio por paso | $2.55\times10^{-2}$ s |
| Inicio registrado | 18 Jun 2026, 16:46:55 |
| Fin registrado | 18 Jun 2026, 18:09:34 |

El tiempo de compilacion no fue medido con una herramienta externa como `time make`, por lo que no se reporta como benchmark cuantitativo. Para una comparacion reproducible de rendimiento, una mejora simple seria ejecutar la compilacion y la simulacion con `/usr/bin/time -p` y guardar esos resultados junto con el log de PLUTO.

---

## 3. Reproducción en Python

### 3.1 Implementación

El script `hall_mhd_harris.py` implementa el pipeline Python usado para reproducir y analizar el problema:
1. **Setup**: condicion inicial de Harris sheet con los mismos parametros de PLUTO.
2. **Cantidades derivadas**: flujo reconectado, corriente $J_z$, divergencia de $\mathbf{B}$ y perfiles 1D.
3. **Visualizacion**: comparacion inicial/final, evolucion temporal y perfiles.
4. **Chequeo numerico**: medicion de $\nabla\cdot\mathbf{B}$ y comparacion directa entre la condicion inicial Python y el snapshot $t=0$ de PLUTO.

El solver evolutivo completo Hall-MHD se deja a PLUTO; la parte Python reproduce la condicion inicial y el post-procesamiento con diferencias finitas, lo cual permite validar setup y diagnosticos sin reimplementar todo el modulo Hall de PLUTO.

### 3.2 Estructura del proyecto y función de cada código

El proyecto se organizo para separar cuatro responsabilidades: archivos de entrada de PLUTO, salidas de la corrida, scripts Python de analisis y reporte. Esta separación evita mezclar el solver externo con los productos generados y permite repetir el flujo de trabajo desde la raiz del repositorio.

| Ruta o archivo | Tipo | Que es | Que hace dentro del proyecto |
|----------------|------|--------|------------------------------|
| `Current_Sheet/init.c` | Código C para PLUTO | Definición física del problema. | Implementa la condición inicial de Harris: $\rho=0.2+\mathrm{sech}^2(y/l)$, $B_x=B_0\tanh(y/l)$, $v_x=v_y=v_z=0$ y una perturbación magnética proporcional a $\Psi_0$. También deja vacías `Analysis()` y `UserDefBoundary()` porque el análisis se hace después en Python. |
| `Current_Sheet/definitions_01.h` | Configuración de compilación PLUTO | Caso principal. | Activa MHD 2D cartesiano, EOS isotérmica, reconstrucción lineal, RK2, `DIV_CLEANING` y `HALL_MHD = EXPLICIT`. Es el archivo que define la física que se compila. |
| `Current_Sheet/definitions_02.h` | Configuración alternativa PLUTO | Variante de comparación. | Se conserva como referencia para una corrida sin Hall o para comparar contra otra configuración. No es el caso principal usado en los resultados finales. |
| `Current_Sheet/pluto_01.ini` | Entrada de ejecución PLUTO | Parámetros numéricos de la corrida principal. | Define dominio, malla, CFL, tiempo final, solver HLL, fronteras, frecuencia de salida y parámetros `ETA`, `WIDTH`, `PSI0`. |
| `Current_Sheet/pluto_02.ini` | Entrada alternativa PLUTO | Parámetros para otra corrida. | Permite repetir pruebas o preparar una comparación futura con otra configuración. |
| `Current_Sheet/pluto_sim/` | Carpeta de corrida | Copia de metadatos y salidas PLUTO. | Guarda `pluto.ini`, `definitions.h`, `init.c`, `vtk.out`, `dbl.out` y, localmente, los dumps `data.*.vtk/.dbl`. Es la fuente que leen los scripts de postproceso. |
| `Current_Sheet/analysis/plot_results.py` | Script Python | Visualización por snapshot. | Lee `pluto_sim/vtk.out`, carga cada snapshot con `pyPLUTO`, calcula presión, $|\mathbf{B}|$ y $J_z$, y genera los paneles `plots/hall_cs_*.png`. |
| `Current_Sheet/analysis/analysis.py` | Script Python | Diagnóstico cuantitativo y comparación. | Lee todos los snapshots, reconstruye la condición inicial analítica, calcula flujo reconectado, $\max |J_z|$, $||\nabla\cdot\mathbf{B}||_2$, errores Python vs PLUTO en $t=0$, CSV y figuras finales del reporte. |
| `Current_Sheet/analysis/figs/` | Salidas de análisis | Figuras y tablas finales. | Contiene `diagnostics.csv`, `initial_condition_errors.csv`, `diagnostics_timeseries.png`, `pluto_final_all_variables.png`, `jz_fieldlines_t55.png` y `python_pluto_initial_comparison.png`. |
| `Current_Sheet/plots/` | Salidas gráficas | Paneles temporales completos. | Contiene una imagen por snapshot desde $t=0$ hasta $t=60$ con todas las variables físicas graficadas. |
| `Current_Sheet/python_reproduction/hall_mhd_harris.py` | Script Python | Reproducción independiente del setup y análisis auxiliar. | Construye la misma malla y condición inicial que PLUTO, define diferencias finitas para $J_z$ y $\nabla\cdot\mathbf{B}$, calcula flujo reconectado y produce gráficas auxiliares. Usa rutas relativas al repositorio. |
| `Current_Sheet/python_reproduction/output/` | Salidas Python | Figuras auxiliares. | Guarda series temporales, comparación inicial/final, evolución de divergencia y perfiles 1D. |
| `Current_Sheet/report/report.md` | Reporte | Documento principal. | Integra teoría, metodología, resultados, comparación, limitaciones, conclusiones y referencias. |
| `Whistler_Waves/` | Configuraciones PLUTO auxiliares | Otro test Hall MHD. | No es el problema elegido para el reporte; se conserva como referencia porque muestra la propagación de ondas whistler, relacionada con la física Hall. |
| `PLUTO/` | Dependencia externa local | Código fuente de PLUTO. | Motor numérico usado para compilar y ejecutar la simulación. No se modifica como parte del análisis del proyecto. |

La simulacion principal se ejecuto con PLUTO hasta $t=60$. Las salidas registradas en `vtk.out` contienen 13 snapshots: $t=0,5,10,\ldots,60$. Los archivos `dbl.out`, `vtk.out`, `pluto.ini`, `definitions.h` e `init.c` dentro de `pluto_sim/` documentan la configuracion exacta usada al momento de correr.

**Lectura rápida de los scripts Python.**

| Script | Entrada | Operaciones principales | Salida |
|--------|---------|-------------------------|--------|
| `analysis/plot_results.py` | Snapshots VTK en `pluto_sim/` | Carga datos con `pyPLUTO`, calcula variables derivadas y grafica cada tiempo. | `plots/hall_cs_*.png` |
| `analysis/analysis.py` | Snapshots VTK y condición inicial analítica | Calcula diagnósticos, errores relativos L2 y figuras comparativas. | `analysis/figs/*.png`, `analysis/figs/*.csv` |
| `python_reproduction/hall_mhd_harris.py` | Parámetros del problema y snapshots PLUTO | Reconstruye setup en Python, calcula $J_z$, flujo reconectado, divergencia y perfiles. | `python_reproduction/output/*.png` |

**Que se corrio.** Primero se configuro y ejecuto PLUTO para la corrida Hall MHD principal. Luego se corrio `analysis/plot_results.py` para producir las figuras por snapshot. Despues se ejecuto `analysis/analysis.py` para calcular los diagnosticos cuantitativos usados en el reporte. Finalmente, `hall_mhd_harris.py` se uso como reproduccion Python del setup y como apoyo para comparar la condicion inicial y las cantidades derivadas.

Los comandos reproducibles usados para reconstruir el flujo de trabajo son:

```bash
# Desde la raiz del repositorio
REPO=$PWD
export PLUTO_DIR="$REPO/PLUTO"   # o la ruta a otra instalacion local de PLUTO

# Configuracion, compilacion y ejecucion con PLUTO
RUN_DIR="$PLUTO_DIR/Test_Problems/MHD/Hall_MHD/Current_Sheet"
cd "$RUN_DIR"
cp "$REPO/Current_Sheet/definitions_01.h" definitions.h
cp "$REPO/Current_Sheet/pluto_01.ini" pluto.ini
cp "$REPO/Current_Sheet/init.c" init.c
python "$PLUTO_DIR/setup.py"
make
./pluto -i pluto.ini | tee pluto_run.log

# Copia de configuracion y bitacoras para reproducibilidad
mkdir -p "$REPO/Current_Sheet/pluto_sim"
cp definitions.h pluto.ini init.c vtk.out dbl.out pluto_run.log "$REPO/Current_Sheet/pluto_sim/"

# Opcional: copiar tambien los dumps crudos si se quiere postprocesar sin rerun
cp data.*.vtk data.*.dbl grid.out "$REPO/Current_Sheet/pluto_sim/"

# Postproceso con pyPLUTO: paneles por snapshot
cd "$REPO"
python Current_Sheet/analysis/plot_results.py

# Diagnosticos cuantitativos y figuras finales del reporte
python Current_Sheet/analysis/analysis.py

# Reproduccion Python del setup y graficas auxiliares
cd Current_Sheet/python_reproduction
python hall_mhd_harris.py
```

En la copia versionada, `pluto_sim/` conserva los archivos minimos para documentar la corrida (`pluto.ini`, `definitions.h`, `init.c`, `vtk.out`, `dbl.out`). Los dumps crudos `data.*.vtk`, `data.*.dbl`, `grid.out`, el ejecutable `pluto` y `pluto_run.log` pueden existir localmente para postproceso, pero son artefactos regenerables y se ignoran en git.

**Que se grafico.** Se graficaron las variables fisicas $\rho$, $P$, $v_x$, $v_y$, $B_x$, $B_y$, $|B|$ y $J_z$ en varios tiempos. Ademas se hicieron curvas temporales del flujo reconectado, de $\max |J_z|$ y del error $||\nabla\cdot B||_2$, una figura especifica de firmas Hall ($B_z$, $(\mathbf{J}\times\mathbf{B})_z/\rho$ y $|\mathbf{v}_e-\mathbf{v}|$), y una figura de corriente con O/X-points anotados. Para la validacion Python vs PLUTO se compararon mapas de $\rho$, $B_x$ y $B_y$ en $t=0$ y se calcularon errores relativos L2.

### 3.3 Resultados del Análisis

| Métrica | Valor |
|---------|-------|
| Flujo reconectado final, $\int_0^{L_x/2}|B_y(x,0)|dx$ | 4.5525 |
| Flujo reconectado firmado final | 2.7213 |
| $|\mathbf{B}|_{max}$ final | 1.2749 |
| $\max(|J_z|)$ inicial | 1.9573 |
| $\max(|J_z|)$ maximo temporal | 3.0819 en $t\approx25$ |
| $\max(|J_z|)$ final | 1.0786 |
| $\nabla\cdot\mathbf{B}$ L2 final | $4.01\times10^{-4}$ |
| $\max |B_z|$ final | 0.2228 |
| $\max |v_z|$ final | 0.4451 |
| $\max |(\mathbf{J}\times\mathbf{B})_z/\rho|$ final | 0.2084 |
| Estimacion local $\psi_{A_z}=|A_z(O)-A_z(X)|$ | 0.9154 |

---

## 4. Analisis Comparativo

### 4.1 Validacion Python vs PLUTO en la condicion inicial

La reproduccion Python de la condicion inicial fue comparada punto a punto contra el snapshot PLUTO en $t=0$. Los errores relativos L2 son:

| Variable | Error relativo L2 | Error absoluto maximo |
|----------|------------------:|----------------------:|
| $\rho$ | $2.47\times10^{-8}$ | $4.67\times10^{-8}$ |
| $v_x$ | 0 | 0 |
| $v_y$ | 0 | 0 |
| $v_z$ | 0 | 0 |
| $B_x$ | $2.54\times10^{-8}$ | $5.96\times10^{-8}$ |
| $B_y$ | $2.60\times10^{-8}$ | $2.33\times10^{-10}$ |
| $B_z$ | 0 | 0 |

Esto confirma que el setup independiente en Python reproduce los campos iniciales usados por PLUTO hasta precision de salida simple. Esta validacion no prueba una evolucion temporal Python, porque Python no actua aqui como solver evolutivo Hall-MHD. Lo que valida es que la condicion inicial, la malla y los diagnosticos derivados fueron reproducidos correctamente. La figura `analysis/figs/python_pluto_initial_comparison.png` muestra los mapas Python, PLUTO y la diferencia para $\rho$, $B_x$ y $B_y$.

### 4.2 Flujo reconectado y corriente

El flujo reconectado crece de $0.040$ en $t=0$ a $4.552$ en $t=60$. El crecimiento se acelera despues de $t\approx20$, alcanza valores cercanos a $4.8$ alrededor de $t=45$, y luego oscila levemente. La corriente maxima aumenta al inicio de la reconexion y alcanza $\max |J_z|\approx3.08$ en $t\approx25$, consistente con la formacion de una capa de corriente intensa.

Una medida mas local de reconexion puede obtenerse mediante el potencial vectorial $A_z$:

$$
\psi(t)=A_z(X)-A_z(O),
$$

donde $X$ y $O$ representan, respectivamente, el punto X de reconexion y el centro de una isla magnetica. Otra opcion fisica es medir la tasa de reconexion con el campo electrico fuera del plano en el punto X:

$$
E_z = -v_xB_y + v_yB_x + \eta J_z + E_{z,\mathrm{Hall}}.
$$

En esta entrega se usa el flujo reconectado integrado sobre el eje medio como diagnostico global porque se calcula directamente desde las salidas VTK. Ademas, se reconstruyo $A_z$ como chequeo local y se obtuvo $\psi_{A_z}\approx0.915$ en $t=60$. El calculo de $E_z$ local queda como extension natural para comparar cuantitativamente Hall MHD contra una corrida ideal o resistiva.

### 4.3 Perfiles 1D y estructura final

Los cortes en $x=0$ muestran la evolucion de las variables a traves de la lamina:
- $B_x$ parte de la configuracion $\tanh(y/l)$ y luego desarrolla estructura alrededor de la zona de reconexion.
- $B_y$ mide el componente reconectado y aumenta despues del onset.
- La densidad se redistribuye por los flujos de salida del evento de reconexion.

### 4.4 Divergencia de B

El error $||\nabla\cdot\mathbf{B}||_2$ inicia en $1.89\times10^{-7}$ y termina en $4.01\times10^{-4}$. Tiene un maximo de $2.24\times10^{-3}$ cerca de $t=30$, todavia pequeno frente a las escalas de campo del problema, y luego decrece. Esto indica que el esquema de divergence cleaning mantiene controlado el error solenoidal durante la corrida.

### 4.5 Comparacion Hall vs no Hall pendiente

La comparacion mas directa para aislar el efecto Hall seria repetir exactamente la misma corrida con `HALL_MHD = NO`, manteniendo dominio, malla, condicion inicial, CFL, solver, tiempo final y frecuencia de salida. Con esa segunda corrida se podrian comparar estas cantidades:

| Cantidad | Que mostraria |
|----------|---------------|
| Flujo reconectado $\int_0^{L_x/2}|B_y(x,0)|dx$ | Si la reconexion crece antes o mas rapido al activar Hall. |
| $\max |J_z|$ | Si la capa de corriente se intensifica o se relaja de forma distinta. |
| Tiempo de onset | El tiempo aproximado en que el flujo reconectado empieza a crecer rapidamente. |
| $\psi_{A_z}=|A_z(O)-A_z(X)|$ | Una comparacion mas fisica de flujo reconectado local. |
| $B_z$ y $(\mathbf{J}\times\mathbf{B})_z/\rho$ | Huellas especificas de la dinamica Hall fuera del plano. |

Esta comparacion no se incluye como resultado cuantitativo porque la corrida no Hall no fue ejecutada todavia. Por lo tanto, el presente trabajo demuestra una evolucion compatible con reconexion bajo Hall MHD y muestra firmas Hall dentro de esa corrida, pero no afirma de forma concluyente que Hall acelere la reconexion frente a MHD ideal o resistiva.

---

## 5. Limitaciones y Posibles Mejoras

### 5.1 Simulación PLUTO

| Aspecto | Limitación | Mejora Posible |
|---------|------------|----------------|
| **Resolución** | 256×128 puede ser insuficiente para capturar la física Hall | Usar mayor resolución (512×256) con AMR de Chombo |
| **Tiempo de cómputo** | ~66 min en CPU para t=60 | Paralelizar con MPI o usar GPU |
| **Hall simplificado** | `ne = rho` (constante) | Implementar $n_e$ físico: $n_e = \rho / (\mu m_p)$ |
| **Solver** | Solo HLL es compatible con Hall | Probar HLLEM para menos difusión numérica |
| **Resistividad** | No se incluyó resistividad explícita | Estudiar el efecto combinado Hall + resistividad |
| **Analisis in situ** | `Analysis()` vacio en init.c | Implementar calculo de flujo reconectado durante la simulacion |
| **Frecuencia de output** | VTK cada 5 unidades | Guardar con mayor frecuencia cerca del onset de reconexión |
| **Configuraciones** | Solo se corrio Hall activo | Correr una variante ideal con `HALL_MHD = NO` para comparacion directa |

### 5.2 Reproducción Python

| Aspecto | Limitación | Mejora Posible |
|---------|------------|----------------|
| **Solver MHD** | Python reproduce setup y analisis, no evoluciona Hall-MHD completo | Implementar HLL/Rusanov 2D con constrained transport |
| **Término Hall** | No implementado en Python | Añadir $-\nabla \times (\mathbf{J} \times \mathbf{B} / n_e e)$ |
| **Upwinding** | Diferencias centradas explotan | Usar esquemas upwind (MUSCL-Hancock, PPM) |
| **Divergencia de B** | Projection method simple | Usar constrained transport (CT) tipo Yee |
| **Rendimiento** | Python puro es lento | Numba/Cython/Mexwell para aceleración |
| **AMR** | Sin adaptación de malla | Implementar refinamiento simple |
| **Validación** | Solo comparación visual | Tests de convergencia, orden del esquema |
| **Documentación** | Comentarios mínimos | Docstrings, referencias a literatura |

### 5.3 Análisis

| Aspecto | Limitación | Mejora Posible |
|---------|------------|----------------|
| **Flujo reconectado** | Integral aproximada | Definición más precisa usando el punto X |
| **Errores L1/L2** | Calculados para la condicion inicial; falta comparacion evolutiva contra otro solver | Comparar Hall vs ideal y hacer estudio de convergencia |
| **Convergencia** | Una sola resolución | Ejecutar mallas 64², 128², 256², 512² |
| **Espectro** | Sin análisis de Fourier | Transformada para ondas whistler |
| **Visualización** | 2D estático | Animaciones, plots interactivos |

### 5.4 Reporte y Documentación

| Aspecto | Mejora |
|---------|--------|
| **Estado del arte** | Revisión más exhaustiva (10-15 referencias) |
| **Literatura Hall** | Birn et al. 2001 (GEM), Huba 2003, Lesur et al. 2014 |
| **Figuras** | Más anotaciones, paneles unificados |
| **Código** | Modularización, tests unitarios, CI |

---

## 6. Conclusiones

1. **La corrida en PLUTO reproduce una evolucion compatible con reconexion magnetica** en una lamina de Harris bajo el modelo Hall MHD, mostrando el flujo reconectado creciendo hasta 4.55 unidades en $t=60$ bajo la definicion $\int_0^{L_x/2}|B_y(x,0)|dx$.

2. **La evolucion observada muestra firmas compatibles con la fisica Hall**: se genera $B_z$ desde una condicion inicial con $B_z=0$, aparece un termino Hall fuera del plano localizado y se observa desacoplamiento efectivo entre velocidad ionica y electronica. Una afirmacion cuantitativa sobre aceleracion frente a MHD ideal requiere correr la variante `HALL_MHD = NO`.

3. **La implementacion en Python** reproduce la condicion inicial con errores relativos L2 de orden $10^{-8}$ frente a PLUTO, calcula cantidades derivadas (corriente, flujo reconectado, divergencia, $A_z$ y diagnosticos Hall) y genera visualizaciones comparativas.

4. **La principal limitación** es doble: falta una corrida no Hall para aislar cuantitativamente el efecto Hall, y Python no actua como solver evolutivo completo. Un solver MHD+Hall desde cero requiere esquemas numéricos avanzados (Riemann, CT) que van más allá del alcance de este proyecto.

5. **Mejoras futuras**: ejecutar la corrida `HALL_MHD = NO`, comparar onset y tasas de reconexion, implementar un calculo local de $E_z$ en el punto X, hacer estudio de convergencia y generar animaciones de la evolucion.

---

## Referencias

- Birn, J., et al. (2001). "Geospace Environmental Modeling (GEM) magnetic reconnection challenge." *JGR*, 106, 3715.
- Harris, E. G. (1962). "On a plasma sheath separating regions of oppositely directed magnetic field." *Il Nuovo Cimento*, 23, 115.
- Huba, J. D. (2003). "Hall Magnetohydrodynamics - A Tutorial." *Space Plasma Simulation*, 615, 166.
- Lesur, G., Kunz, M. W., & Fromang, S. (2014). "Thanatology in protoplanetary discs." *A&A*, 566, A56.
- Mignone, A., et al. (2012). "The PLUTO Code for Adaptive Mesh Computations in Astrophysical Fluid Dynamics." *ApJS*, 198, 7.
- Viganò, D., Pons, J. A., & Miralles, J. A. (2012). "A new code for the Hall-driven magnetic evolution of neutron stars." *CPC*, 183, 2042.
