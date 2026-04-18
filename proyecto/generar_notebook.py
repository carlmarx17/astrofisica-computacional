"""
generar_notebook.py  —  Construye el cuaderno Jupyter de forma programática.

Escribimos el cuaderno como un script de Python en lugar de editar el JSON
del .ipynb directamente, porque así el código fuente es legible y los
diff de git son realmente útiles.

Ejecuta esto cada vez que cambies el contenido del cuaderno:
    python generar_notebook.py
"""

import nbformat as nbf

cuaderno = nbf.v4.new_notebook()


# ─────────────────────────────────────────────────────────────────────────────
# Portada
# ─────────────────────────────────────────────────────────────────────────────

portada = r"""# 🌞 Viento Solar de Parker — Simulación Numérica
*Astrofísica Computacional · Proyecto 1*

---

Modelamos cómo la corona del Sol sopla plasma hacia afuera a cientos de km/s —
y por qué debe volverse supersónico antes de llegar a la Tierra.

El actor principal es la ecuación de momento de Parker de 1958 (con el término de advección):

$$
v \frac{dv}{dr} = \frac{2 v_c^2}{r} - \frac{G M_\odot}{r^2}
$$

donde $v_c = \sqrt{k_B T / \mu m_p}$ es la velocidad del sonido isotérmica y
$r_c = G M_\odot / 2v_c^2$ es el **radio crítico** — el punto de no retorno.

Implementamos dos integradores manuales (Euler y RK4) y los validamos
contra la solución semi-analítica obtenida por búsqueda de raíces.
"""


# ─────────────────────────────────────────────────────────────────────────────
# Importaciones
# ─────────────────────────────────────────────────────────────────────────────

codigo_importaciones = r"""%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12, 'figure.dpi': 110})

# Nuestros módulos personalizados (todos en el mismo directorio)
from constantes import *
from analitico import donde_se_vuelve_sonico, radio_sonico, cazar_velocidad_en_todo_el_espacio
from solucionadores_numericos import lanzar_viento_solar
"""


# ─────────────────────────────────────────────────────────────────────────────
# Parte 1 — Solución Analítica
# ─────────────────────────────────────────────────────────────────────────────

md_parte1 = r"""## Parte 1 · La Referencia Analítica

Empezamos resolviendo el problema *exactamente* (bueno, semi-analíticamente con búsqueda de raíces).
La ley de conservación integrada de Parker nos da una ecuación implícita para $v(r)$:

$$
\left(\frac{v}{v_c}\right)^2 - \ln\left(\frac{v}{v_c}\right)^2
= 4\ln\frac{r}{r_c} + \frac{4 r_c}{r} - 3
$$

Llamamos a `cazar_velocidad_en_todo_el_espacio()` que aplica el método de Brent en cada radio para
encontrar la única solución transónica (subsónica para $r < r_c$, supersónica para $r > r_c$).
"""

codigo_parte1 = r"""# Malla de 1000 radios de 1 a 100 radios solares (en metros)
malla_r = np.linspace(1.0, 100.0, 1_000) * RADIO_SOL

# ── Valores del punto crítico ────────────────────────────────────────────────
v_c = donde_se_vuelve_sonico(TEMPERATURA_BASE)     # velocidad crítica  [m/s]
r_c = radio_sonico(TEMPERATURA_BASE)               # radio crítico [m]

print(f"Velocidad crítica  v_c = {v_c/1e3:.1f} km/s   (esperado ~91 km/s)")
print(f"Radio crítico      r_c = {r_c/RADIO_SOL:.2f} R☉    (esperado ~5.8 R☉)")

# ── Cazar la solución transónica en cada punto de la malla ───────────────────
v_exacta = cazar_velocidad_en_todo_el_espacio(malla_r, TEMPERATURA_BASE)

# ── Graficar ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(malla_r / RADIO_SOL, v_exacta / 1e3, 'k-', lw=2.5, label='Analítica (rama transónica)')
ax.plot(r_c / RADIO_SOL, v_c / 1e3, 'ro', ms=10, zorder=5, label='Punto crítico $(r_c, v_c)$')

# Líneas de referencia punteadas en el punto sónico
ax.axhline(v_c / 1e3, color='gray', ls='--', alpha=0.6)
ax.axvline(r_c / RADIO_SOL, color='gray', ls='--', alpha=0.6)

# Etiquetar las dos regiones
ax.text( 2, v_c/1e3 - 15, 'Subsónico',    color='steelblue', fontsize=11)
ax.text(10, v_c/1e3 + 20, 'Supersónico',  color='firebrick', fontsize=11)

ax.set_xlabel(r'Radio  $r$  ($R_\odot$)')
ax.set_ylabel(r'Velocidad  $v$  (km s$^{-1}$)')
ax.set_title(r'Solución Analítica del Viento de Parker  ($T = 10^6$ K)')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()
"""


# ─────────────────────────────────────────────────────────────────────────────
# Parte 2 — Integradores Numéricos
# ─────────────────────────────────────────────────────────────────────────────

md_parte2 = r"""## Parte 2 · Integración Numérica — Pasos de Bebé vs. Los Cuatro Elegantes

Ahora integramos la EDO numéricamente sin usar solucionadores incorporados.

| Método | Orden | Tamaño de paso | Notas |
|--------|-------|----------------|-------|
| Euler (`pasos_de_bebe_euler`) | 1º | ~10⁵ m | Simple pero se desvía cerca de regiones rígidas |
| RK4 (`los_cuatro_elegantes_rk4`) | 4º | ~10⁶ m | Mucho más preciso con pasos 10× mayores |

Ambos métodos parten de una pequeña perturbación $\varepsilon = 10^{-3}$ lejos de
$(r_c, v_c)$ para evitar la singularidad $0/0$, luego integran hacia afuera (supersónico)
y hacia adentro (subsónico) por separado.
"""

codigo_parte2 = r"""# ── Resolver con ambos métodos ──────────────────────────────────────────────
r_rk4, v_rk4, rho_rk4 = lanzar_viento_solar(TEMPERATURA_BASE, metodo='rk4',   tamano_paso=1e6)
r_eu,  v_eu,  rho_eu  = lanzar_viento_solar(TEMPERATURA_BASE, metodo='euler',  tamano_paso=1e5)

# ── Graficar la comparación ───────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(malla_r / RADIO_SOL, v_exacta / 1e3, 'k-',  lw=3,   label='Analítica exacta')
ax.plot(r_rk4,  v_rk4,                       'b--', lw=2,   label='RK4   (h = 1 000 km)')
ax.plot(r_eu,   v_eu,                         'r:',  lw=2,   label='Euler (h =   100 km)')

ax.set_xlabel(r'Radio  $r$  ($R_\odot$)')
ax.set_ylabel(r'Velocidad  $v$  (km s$^{-1}$)')
ax.set_title('Euler vs RK4 — ¿Qué tan cerca llegamos a la respuesta exacta?')
ax.set_xlim(1, 100)
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# ── Errores relativos en r = 50 R☉ ────────────────────────────────────────────
v_exacta_50 = np.interp(50.0, malla_r / RADIO_SOL, v_exacta / 1e3)
v_rk4_50    = np.interp(50.0, r_rk4, v_rk4)
v_eu_50     = np.interp(50.0, r_eu,  v_eu)

print(f"Velocidad exacta en 50 R☉    = {v_exacta_50:.4f} km/s")
print(f"Error relativo de RK4         = {abs(v_rk4_50 - v_exacta_50)/v_exacta_50:.2e}")
print(f"Error relativo de Euler       = {abs(v_eu_50  - v_exacta_50)/v_exacta_50:.2e}")
"""


# ─────────────────────────────────────────────────────────────────────────────
# Parte 3 — Sensibilidad a la Temperatura
# ─────────────────────────────────────────────────────────────────────────────

md_parte3 = r"""## Parte 3 · Sensibilidad a la Temperatura — ¿Qué pasa si la corona está más caliente?

Mayor temperatura → mayor velocidad del sonido → el punto sónico se acerca al Sol
→ el viento es más rápido y menos denso a una distancia dada.

Probamos $T \in \{0.5,\ 1,\ 2\} \times 10^6$ K para ver qué tan sensibles son
la curva de aceleración y el perfil de densidad a la temperatura coronal.
"""

codigo_parte3 = r"""temperaturas = [0.5e6, 1e6, 2e6]
paleta       = ['forestgreen', 'steelblue', 'darkorchid']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

for T, color in zip(temperaturas, paleta):
    r_sol, v_sol, rho_sol = lanzar_viento_solar(T, metodo='rk4', tamano_paso=1e6)

    # Marcar el punto crítico para cada temperatura
    r_sonico = radio_sonico(T) / RADIO_SOL
    v_sonico = donde_se_vuelve_sonico(T) / 1e3

    etiqueta = f'T = {T/1e6:.1f} MK'
    ax1.plot(r_sol, v_sol, color=color, lw=2,   label=etiqueta)
    ax1.plot(r_sonico, v_sonico, 'o', color=color, ms=7)

    ax2.loglog(r_sol, rho_sol, color=color, lw=2, label=etiqueta)

ax1.set_xlabel(r'Radio  $r$  ($R_\odot$)')
ax1.set_ylabel(r'Velocidad  $v$  (km s$^{-1}$)')
ax1.set_title('Velocidad del Viento a Distintas Temperaturas')
ax1.legend()
ax1.grid(alpha=0.3)

ax2.set_xlabel(r'Radio  $r$  ($R_\odot$)')
ax2.set_ylabel(r'Densidad  $\rho$  (kg m$^{-3}$)')
ax2.set_title('Densidad de Masa (log-log)')
ax2.legend()
ax2.grid(alpha=0.3)

plt.suptitle('Estudio de Sensibilidad a la Temperatura', fontsize=14, y=1.01)
plt.tight_layout()
plt.show()

# ── Condiciones en 1 UA (≈215 R☉) ─────────────────────────────────────────────
print("Condiciones del viento en 1 UA (215 R☉):")
print(f"{'Temperatura':>18}  {'Velocidad (km/s)':>16}")
print("-" * 38)
for T in temperaturas:
    r_sol, v_sol, _ = lanzar_viento_solar(T, r_max=220, metodo='rk4', tamano_paso=2e6)
    v_en_1ua = np.interp(215.0, r_sol, v_sol)
    print(f"{T:.2e} K            {v_en_1ua:>12.1f}")
"""


# ─────────────────────────────────────────────────────────────────────────────
# Parte 4 — Variación de Tasa de Pérdida de Masa
# ─────────────────────────────────────────────────────────────────────────────

md_parte4 = r"""## Parte 4 · Variación de Tasa de Pérdida de Masa — ¿Más o menos viento?

La ecuación de momento no contiene $\dot{M}$, así que el perfil de velocidad
$v(r)$ es **completamente independiente** de la tasa de pérdida de masa.
Solo escala la densidad: $\rho \propto \dot{M}$.

Verificamos esto corriendo con $\dot{M}/2$, $\dot{M}$, y $2\dot{M}$.
"""

codigo_parte4 = r"""multiplicadores_flujo = [0.5, 1.0, 2.0]
paleta                = ['darkorange', 'black', 'mediumblue']

# Una sola corrida de velocidad es suficiente — v no depende de Ṁ
r_base, v_base, _ = lanzar_viento_solar(TEMPERATURA_BASE, r_max=50, metodo='rk4', tamano_paso=1e6)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Izquierda: velocidad (sin cambios con Ṁ)
ax1.plot(r_base, v_base, 'k-', lw=2.5)
ax1.set_xlabel(r'Radio  $r$  ($R_\odot$)')
ax1.set_ylabel(r'Velocidad  $v$  (km s$^{-1}$)')
ax1.set_title(r'Velocidad — independiente de $\dot{M}$')
ax1.set_xlim(1, 50)
ax1.grid(alpha=0.3)
ax1.text(5, max(v_base)*0.5, "La misma curva para TODOS los valores de Ṁ",
         fontsize=10, color='gray', style='italic')

# Derecha: densidad (escala linealmente con Ṁ)
for f, color in zip(multiplicadores_flujo, paleta):
    # Recalcular densidad con el Ṁ escalado
    rho_escalada = (f * FLUJO_MASA_BASE) / (
        4.0 * np.pi * (r_base * RADIO_SOL)**2 * (v_base * 1e3)
    )
    ax2.loglog(r_base, rho_escalada, color=color, lw=2, label=rf'$\dot{{M}}={f}\times\dot{{M}}_0$')

ax2.set_xlabel(r'Radio  $r$  ($R_\odot$)')
ax2.set_ylabel(r'Densidad  $\rho$  (kg m$^{-3}$)')
ax2.set_title(r'La densidad escala linealmente con $\dot{M}$')
ax2.legend()
ax2.set_xlim(1, 50)
ax2.grid(alpha=0.3)

plt.suptitle('Estudio de Tasa de Pérdida de Masa', fontsize=14, y=1.01)
plt.tight_layout()
plt.show()
"""


# ─────────────────────────────────────────────────────────────────────────────
# Parte 5 — Extensión con Calentamiento
# ─────────────────────────────────────────────────────────────────────────────

md_parte5 = r"""## Parte 5 · Extensión — ¿Podemos empujar el viento más rápido con calor?

Los modelos reales del viento solar incluyen calentamiento coronal más allá del isotérmico.
Agregamos una fuente de calor exponencial simple:

$$Q(r) = Q_0 \, e^{-r/H}, \quad H = 5\,R_\odot$$

Este término extra $Q(r)/v$ se agrega al numerador de $dv/dr$.
Ajustamos $Q_0$ para lograr aproximadamente **+20% de velocidad terminal** respecto a la base.
"""

codigo_parte5 = r"""# ── Definir la función de calentamiento ───────────────────────────────────────
escala_H = 5.0 * RADIO_SOL   # altura característica en metros

# Q_0 ajustado para producir ~20% de velocidad terminal extra
Q0 = 0.5 * (donde_se_vuelve_sonico(TEMPERATURA_BASE)**2 / RADIO_SOL)

def bomba_de_calor_coronal(r):
    # Explosión exponencial de energía depositada cerca del Sol
    return Q0 * np.exp(-r / escala_H)

# ── Correr ambas versiones ────────────────────────────────────────────────────
r_caliente, v_caliente, _ = lanzar_viento_solar(
    TEMPERATURA_BASE, metodo='rk4', tamano_paso=1e6,
    funcion_calor=bomba_de_calor_coronal
)

# ── Comparar ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(r_rk4,      v_rk4,      color='steelblue', lw=2.5, label='Isotérmico base')
ax.plot(r_caliente, v_caliente, color='firebrick',  lw=2.5, ls='--',
        label='Con calentamiento exponencial Q(r)')

ax.set_xlabel(r'Radio  $r$  ($R_\odot$)')
ax.set_ylabel(r'Velocidad  $v$  (km s$^{-1}$)')
ax.set_title('Efecto del Calentamiento Coronal en la Velocidad Terminal')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# ── Imprimir el aumento de velocidad ──────────────────────────────────────────
v_base_100     = np.interp(100.0, r_rk4,      v_rk4)
v_caliente_100 = np.interp(100.0, r_caliente, v_caliente)
aumento        = (v_caliente_100 - v_base_100) / v_base_100 * 100.0

print(f"Velocidad terminal — base          : {v_base_100:.1f} km/s")
print(f"Velocidad terminal — con Q(r)      : {v_caliente_100:.1f} km/s")
print(f"Aumento de velocidad               : {aumento:+.1f} %")
"""


# ─────────────────────────────────────────────────────────────────────────────
# Armar y escribir el cuaderno
# ─────────────────────────────────────────────────────────────────────────────

cuaderno['cells'] = [
    nbf.v4.new_markdown_cell(portada),
    nbf.v4.new_code_cell(codigo_importaciones),
    nbf.v4.new_markdown_cell(md_parte1),
    nbf.v4.new_code_cell(codigo_parte1),
    nbf.v4.new_markdown_cell(md_parte2),
    nbf.v4.new_code_cell(codigo_parte2),
    nbf.v4.new_markdown_cell(md_parte3),
    nbf.v4.new_code_cell(codigo_parte3),
    nbf.v4.new_markdown_cell(md_parte4),
    nbf.v4.new_code_cell(codigo_parte4),
    nbf.v4.new_markdown_cell(md_parte5),
    nbf.v4.new_code_cell(codigo_parte5),
]

with open('viento_solar_parker.ipynb', 'w') as archivo:
    nbf.write(cuaderno, archivo)

print("✅  viento_solar_parker.ipynb escrito exitosamente!")
