# =============================================================================
# solucionadores_numericos.py  —  "Haciendo la integral a mano, como Parker quería"
# =============================================================================
# Implementamos DOS integradores manuales de EDO para la ecuación del viento de Parker:
#
#   1.  Pasos de bebé (Euler)      — primer orden, simple pero necesita pasos pequeños
#   2.  Los cuatro elegantes (RK4) — cuarto orden, preciso con pasos más grandes
#
# Ambos métodos deben manejar una singularidad molesta en el punto crítico r_c
# donde dv/dr = 0/0.  Nuestro truco: no empezar ahí.  En cambio, nos alejamos
# ligeramente de r_c (usando un ε pequeño) e integramos
# hacia afuera (rama supersónica) y hacia adentro (rama subsónica) por separado,
# luego unimos las dos piezas en una solución continua.
# =============================================================================

import numpy as np

from constantes import MASA_SOL, GRAVEDAD, RADIO_SOL, FLUJO_MASA_BASE


# ─────────────────────────────────────────────────────────────────────────────
# La EDO en sí: dv/dr
# ─────────────────────────────────────────────────────────────────────────────

def aceleracion_del_viento(r, v, velocidad_sonica, calentamiento=0.0):
    """
    Calcula el gradiente de velocidad dv/dr en la posición r.

    Esto proviene de combinar la ecuación de momento
        ρ v dv/dr = -dP/dr - ρ G M / r²
    con la ecuación de estado isotérmica P = (k_B T / µ m_p) ρ.

    Después de sustituir y simplificar:

        dv/dr = [ 2 v_c² / r  −  G M / r²  +  Q(r) ]
                ─────────────────────────────────────
                        v  −  v_c² / v

    Nota: El denominador se anula en v = v_c (punto sónico), por eso
    no podemos iniciar la integración ahí.

    Parámetros
    ----------
    r              : float  — posición radial [m]
    v              : float  — velocidad del viento en r [m/s]
    velocidad_sonica: float — velocidad crítica v_c [m/s]
    calentamiento  : float  — calentamiento extra opcional Q(r) [m/s²], por defecto 0

    Regresa
    -------
    float
        dv/dr en el punto (r, v) dado.
    """
    numerador   = (2.0 * velocidad_sonica**2 / r) - (GRAVEDAD * MASA_SOL / r**2) + calentamiento
    denominador = v - (velocidad_sonica**2 / v)
    return numerador / denominador


# ─────────────────────────────────────────────────────────────────────────────
# Integrador 1: Pasos de bebé (Euler explícito)
# ─────────────────────────────────────────────────────────────────────────────

def pasos_de_bebe_euler(r, v, dr, velocidad_sonica, calor_en_r=0.0):
    """
    Avanza la solución un paso usando el método de Euler explícito.

    Literalmente: v_nueva = v + dv/dr · Δr
    Es simple, intuitivo, y tiende a desviarse o explotar cerca de regiones
    rígidas a menos que uses Δr muy pequeño.  Piénsalo como "suma de Riemann
    por la izquierda en el tiempo."

    Parámetros
    ----------
    r              : float — posición radial actual [m]
    v              : float — velocidad actual del viento [m/s]
    dr             : float — tamaño del paso [m] (puede ser negativo para integrar hacia adentro)
    velocidad_sonica: float — v_c [m/s]
    calor_en_r     : float — Q(r) en el punto actual [m/s²]

    Regresa
    -------
    float
        Velocidad actualizada en r + dr [m/s].
    """
    pendiente = aceleracion_del_viento(r, v, velocidad_sonica, calor_en_r)
    return v + pendiente * dr


# ─────────────────────────────────────────────────────────────────────────────
# Integrador 2: Los cuatro elegantes (RK4)
# ─────────────────────────────────────────────────────────────────────────────

def los_cuatro_elegantes_rk4(r, v, dr, velocidad_sonica, funcion_calor=None):
    """
    Avanza la solución un paso usando Runge-Kutta de 4to orden.

    RK4 evalúa la pendiente en cuatro puntos estratégicos dentro de cada paso
    y toma un promedio ponderado — es mucho más preciso que Euler
    con el mismo (o mayor) tamaño de paso.  Piénsalo como "regla de Simpson
    en el tiempo."

    Las cuatro pendientes:
        k1 = f(r,          v          )   — pendiente al inicio
        k2 = f(r + dr/2,   v + k1·dr/2)   — pendiente a la mitad (usando k1)
        k3 = f(r + dr/2,   v + k2·dr/2)   — pendiente a la mitad (usando k2)
        k4 = f(r + dr,     v + k3·dr  )   — pendiente al final

    Actualización final: v_nueva = v + (dr/6) · (k1 + 2k2 + 2k3 + k4)

    Parámetros
    ----------
    r             : float    — posición radial actual [m]
    v             : float    — velocidad actual del viento [m/s]
    dr            : float    — tamaño del paso [m]
    velocidad_sonica: float  — v_c [m/s]
    funcion_calor : callable — Q(r) opcional; se llama como funcion_calor(r)

    Regresa
    -------
    float
        Velocidad actualizada en r + dr [m/s].
    """
    # Función auxiliar para evaluar el calentamiento opcional de forma segura
    def Q(radio):
        return funcion_calor(radio) if funcion_calor else 0.0

    k1 = aceleracion_del_viento(r,        v,              velocidad_sonica, Q(r))
    k2 = aceleracion_del_viento(r + dr/2, v + k1*dr/2,   velocidad_sonica, Q(r + dr/2))
    k3 = aceleracion_del_viento(r + dr/2, v + k2*dr/2,   velocidad_sonica, Q(r + dr/2))
    k4 = aceleracion_del_viento(r + dr,   v + k3*dr,     velocidad_sonica, Q(r + dr))

    return v + (dr / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)


# ─────────────────────────────────────────────────────────────────────────────
# Solucionador principal: ¡lanzar el viento solar!
# ─────────────────────────────────────────────────────────────────────────────

def lanzar_viento_solar(
    temperatura,
    r_min=1.5,
    r_max=100.0,
    tamano_paso=1e6,
    metodo='rk4',
    funcion_calor=None
):
    """
    Resuelve la EDO del viento solar de Parker y devuelve v(r) y ρ(r).

    Estrategia
    ----------
    Como dv/dr explota exactamente en el punto crítico (r_c, v_c),
    no podemos simplemente empezar ahí.  En cambio:

      1. Empezamos en r_c + ε  (justo supersónico) y marchamos hacia r_max.
      2. Empezamos en r_c - ε  (justo subsónico)   y marchamos hacia r_min.
      3. Unimos ambas ramas en r_c.

    Parámetros
    ----------
    temperatura   : float    — temperatura de la corona [K]
    r_min         : float    — frontera interna [radios solares], por defecto 1.5 R☉
    r_max         : float    — frontera externa [radios solares], por defecto 100 R☉
    tamano_paso   : float    — tamaño del paso de integración [m], por defecto 1e6 m
    metodo        : str      — 'rk4' (recomendado) o 'euler' (educativo)
    funcion_calor : callable — calentamiento coronal opcional Q(r) [m/s²]

    Regresa
    -------
    r_en_radios_solares : np.ndarray — radios en radios solares  (R☉)
    v_en_km_por_s       : np.ndarray — velocidad del viento en km/s
    densidad            : np.ndarray — densidad de masa en kg/m³
    """
    from analitico import donde_se_vuelve_sonico, radio_sonico

    v_c = donde_se_vuelve_sonico(temperatura)
    r_c = radio_sonico(temperatura)

    # Convertir fronteras de radios solares a metros
    r_min_m = r_min * RADIO_SOL
    r_max_m = r_max * RADIO_SOL

    # Un pequeño empujón para no empezar exactamente en la singularidad
    epsilon = 1e-3

    # ── Rama A: Marchar hacia afuera (supersónica) ────────────────────────────

    r_afuera, v_afuera = [], []
    r_ahora = r_c * (1.0 + epsilon)   # justo pasado el punto sónico
    v_ahora = v_c * (1.0 + epsilon)   # justo encima de la velocidad crítica
    dr = +tamano_paso                  # paso positivo (hacia afuera)

    r_afuera.append(r_ahora)
    v_afuera.append(v_ahora)

    while r_ahora < r_max_m:
        # Recortar el último paso para aterrizar exactamente en r_max
        if r_ahora + dr > r_max_m:
            dr = r_max_m - r_ahora

        if metodo == 'rk4':
            v_ahora = los_cuatro_elegantes_rk4(r_ahora, v_ahora, dr, v_c, funcion_calor)
        else:
            calor_aqui = funcion_calor(r_ahora) if funcion_calor else 0.0
            v_ahora = pasos_de_bebe_euler(r_ahora, v_ahora, dr, v_c, calor_aqui)

        r_ahora += dr
        r_afuera.append(r_ahora)
        v_afuera.append(v_ahora)

    # ── Rama B: Marchar hacia adentro (subsónica) ─────────────────────────────

    r_adentro, v_adentro = [], []
    r_ahora = r_c * (1.0 - epsilon)   # justo dentro del punto sónico
    v_ahora = v_c * (1.0 - epsilon)   # justo debajo de la velocidad crítica
    dr = -tamano_paso                  # paso negativo (hacia adentro)

    r_adentro.append(r_ahora)
    v_adentro.append(v_ahora)

    while r_ahora > r_min_m:
        # Recortar para no pasar de r_min
        if r_ahora + dr < r_min_m:
            dr = r_min_m - r_ahora

        if metodo == 'rk4':
            v_ahora = los_cuatro_elegantes_rk4(r_ahora, v_ahora, dr, v_c, funcion_calor)
        else:
            calor_aqui = funcion_calor(r_ahora) if funcion_calor else 0.0
            v_ahora = pasos_de_bebe_euler(r_ahora, v_ahora, dr, v_c, calor_aqui)

        r_ahora += dr
        r_adentro.append(r_ahora)
        v_adentro.append(v_ahora)

    # ── Unir todo: subsónica ← sónica → supersónica ──────────────────────────

    # La rama interna se construyó de r_c hacia r_min; hay que voltearla
    r_adentro = r_adentro[::-1]
    v_adentro = v_adentro[::-1]

    r_total = np.concatenate((r_adentro, [r_c], r_afuera))
    v_total = np.concatenate((v_adentro, [v_c], v_afuera))

    # ── Densidad de la ecuación de continuidad: ρ = Ṁ / (4π r² v) ────────────

    densidad = FLUJO_MASA_BASE / (4.0 * np.pi * r_total**2 * v_total)

    # ── Convertir a unidades amigables para la salida ─────────────────────────

    r_en_radios_solares = r_total / RADIO_SOL   # metros → radios solares
    v_en_km_por_s       = v_total / 1000.0      # m/s    → km/s

    return r_en_radios_solares, v_en_km_por_s, densidad
