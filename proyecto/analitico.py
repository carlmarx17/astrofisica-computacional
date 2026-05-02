# =============================================================================
# analitico.py  —  "El universo tiene una respuesta exacta, vamos a encontrarla"
# =============================================================================
# Este módulo contiene las soluciones semi-analíticas para el viento solar de Parker.
# La idea clave: Parker derivó una ecuación implícita que relaciona v y r.
# La giramos y preguntamos "¿para cada r, cuál es v?" usando búsqueda de raíces.
#
# Ecuación clave (Parker 1958, integrada una vez):
#
#   (v/v_c)² - ln(v/v_c)² = 4·ln(r/r_c) + 4·(r_c/r) - 3
#
# Donde r_c y v_c son el punto "crítico" (sónico) donde el viento se vuelve supersónico.
# =============================================================================

import numpy as np
from scipy.optimize import brentq

from constantes import MASA_SOL, GRAVEDAD, BOLTZMANN, MASA_PROTON, PESO_MOLECULAR, RADIO_SOL


# ─────────────────────────────────────────────────────────────────────────────
# Paso 1: Encontrar los números mágicos (el punto crítico)
# ─────────────────────────────────────────────────────────────────────────────

def donde_se_vuelve_sonico(temperatura):
    """
    Calcula la velocidad sónica (crítica) v_c.

    Esta es la velocidad a la que el viento solar cruza de subsónico a supersónico.
    Resulta ser la velocidad del sonido isotérmica a la temperatura dada:

        v_c = sqrt( k_B · T / (µ · m_p) )

    Parámetros
    ----------
    temperatura : float
        Temperatura de la corona en Kelvin.

    Regresa
    -------
    float
        Velocidad crítica v_c en m/s. Para T=1 MK es ~91 km/s.
    """
    return np.sqrt(BOLTZMANN * temperatura / (PESO_MOLECULAR * MASA_PROTON))


def radio_sonico(temperatura):
    """
    Encuentra la distancia radial donde el viento solar alcanza la velocidad sónica.

    En este "radio crítico" r_c, el equilibrio entre la gravedad del Sol
    y el gradiente de presión fuerza al viento a ir exactamente v = v_c.

        r_c = G · M_sol / (2 · v_c²)

    Parámetros
    ----------
    temperatura : float
        Temperatura de la corona en Kelvin.

    Regresa
    -------
    float
        Radio crítico r_c en metros. Para T=1 MK es ~5-6 radios solares.
    """
    v_c = donde_se_vuelve_sonico(temperatura)
    return GRAVEDAD * MASA_SOL / (2.0 * v_c**2)


# ─────────────────────────────────────────────────────────────────────────────
# Paso 2: La ecuación trascendental implícita — confía en las matemáticas
# ─────────────────────────────────────────────────────────────────────────────

def ley_de_conservacion_de_parker(v, r, v_c, r_c):
    """
    El residuo de la ecuación integrada del viento de Parker.

    La ecuación del viento de Parker, integrada una vez, da una relación
    similar a una ley de conservación entre v y r. La reacomodamos
    para que la velocidad "correcta" sea la que hace esto igual a cero:

        f(v) = (v/v_c)² - ln(v/v_c)² - 4·ln(r/r_c) - 4·(r_c/r) + 3  ≡ 0

    Parámetros
    ----------
    v   : float  — velocidad de prueba [m/s]
    r   : float  — posición radial [m]
    v_c : float  — velocidad crítica [m/s]
    r_c : float  — radio crítico [m]

    Regresa
    -------
    float
        El residuo. Una solución perfecta devuelve exactamente 0.
    """
    v_norm = v / v_c
    r_norm = r / r_c

    # Lado izquierdo de la ecuación: cuadrático + logaritmo en número de Mach
    lado_izquierdo = v_norm**2 - np.log(v_norm**2)

    # Lado derecho: términos geométricos que dependen de la posición
    lado_derecho = 4.0 * np.log(r_norm) + 4.0 / r_norm - 3.0

    return lado_izquierdo - lado_derecho


def derivada_ley_de_conservacion_de_parker(v, v_c):
    """
    Derivada de la ecuación implícita de Parker respecto a la velocidad.

    La usamos para aplicar Newton-Raphson:

        v_{n+1} = v_n - f(v_n) / f'(v_n)
    """
    return (2.0 * v / v_c**2) - (2.0 / v)


def encontrar_velocidad_por_newton_raphson(r, v_c, r_c, estimado_inicial,
                                           tolerancia=1e-10, max_iter=100):
    """
    Resuelve la ecuación implícita de Parker en un radio fijo con Newton-Raphson.

    Regresa NaN si el método pisa una región no física (v <= 0), si la derivada
    se vuelve demasiado pequeña o si no converge en `max_iter`.
    """
    v = float(estimado_inicial)

    for _ in range(max_iter):
        if v <= 0.0:
            return np.nan

        residuo = ley_de_conservacion_de_parker(v, r, v_c, r_c)
        derivada = derivada_ley_de_conservacion_de_parker(v, v_c)

        if abs(derivada) < 1e-14:
            return np.nan

        paso = residuo / derivada
        v_nueva = v - paso

        if v_nueva <= 0.0:
            return np.nan

        if abs(v_nueva - v) <= tolerancia * max(1.0, abs(v_nueva)):
            return v_nueva

        v = v_nueva

    return np.nan


# ─────────────────────────────────────────────────────────────────────────────
# Paso 3: Cazar la velocidad del viento en cada punto del espacio
# ─────────────────────────────────────────────────────────────────────────────

def cazar_velocidad_en_todo_el_espacio(arreglo_r, temperatura):
    """
    Resuelve v(r) analíticamente en cada radio del arreglo_r.

    Usamos el método de búsqueda de raíces de Brent sobre `ley_de_conservacion_de_parker`.
    La parte difícil: hay DOS soluciones en cada r (la rama subsónica
    y la supersónica). Elegimos la correcta limitando el intervalo de búsqueda:
      - Para r < r_c → rama subsónica, buscamos v < v_c
      - Para r > r_c → rama supersónica, buscamos v > v_c

    Parámetros
    ----------
    arreglo_r   : np.ndarray  — radios en metros
    temperatura : float       — temperatura de la corona en Kelvin

    Regresa
    -------
    np.ndarray
        Velocidad del viento en m/s en cada radio. NaN si la búsqueda falla.
    """
    v_c = donde_se_vuelve_sonico(temperatura)
    r_c = radio_sonico(temperatura)

    velocidades = np.zeros_like(arreglo_r, dtype=float)

    for i, r in enumerate(arreglo_r):

        if r == r_c:
            # Exactamente en el punto crítico — ¡fácil!
            velocidades[i] = v_c

        elif r < r_c:
            # Región subsónica: el plasma sigue arrastrándose (v < v_c)
            try:
                velocidades[i] = brentq(
                    ley_de_conservacion_de_parker,
                    a=1e-10,          # límite inferior: velocidad prácticamente cero
                    b=v_c - 1e-5,     # límite superior: justo debajo de la velocidad sónica
                    args=(r, v_c, r_c)
                )
            except ValueError:
                velocidades[i] = np.nan  # no debería pasar, pero por si acaso

        else:
            # Región supersónica: el plasma sale disparado hacia afuera (v > v_c)
            try:
                velocidades[i] = brentq(
                    ley_de_conservacion_de_parker,
                    a=v_c + 1e-5,     # límite inferior: justo encima de la velocidad sónica
                    b=100.0 * v_c,    # límite superior: generoso físicamente
                    args=(r, v_c, r_c)
                )
            except ValueError:
                velocidades[i] = np.nan

    return velocidades


def cazar_velocidad_en_todo_el_espacio_newton_raphson(arreglo_r, temperatura,
                                                      tolerancia=1e-10, max_iter=100):
    """
    Resuelve v(r) con Newton-Raphson usando la solución previa como semilla.

    Es más simple conceptualmente que Brent, pero también menos robusto:
    depende de una buena semilla inicial y evita el punto crítico usando
    el valor exacto v = v_c cuando r = r_c.
    """
    v_c = donde_se_vuelve_sonico(temperatura)
    r_c = radio_sonico(temperatura)

    arreglo_r = np.asarray(arreglo_r, dtype=float)
    velocidades = np.zeros_like(arreglo_r, dtype=float)

    indices_subsonicos = np.where(arreglo_r < r_c)[0]
    indices_supersonicos = np.where(arreglo_r > r_c)[0]
    indices_criticos = np.where(np.isclose(arreglo_r, r_c, rtol=0.0, atol=1e-12 * r_c))[0]

    velocidades[indices_criticos] = v_c

    semilla_subsonica = 0.95 * v_c
    for i in indices_subsonicos[::-1]:
        velocidad = encontrar_velocidad_por_newton_raphson(
            arreglo_r[i], v_c, r_c, semilla_subsonica,
            tolerancia=tolerancia, max_iter=max_iter
        )
        velocidades[i] = velocidad
        if np.isfinite(velocidad):
            semilla_subsonica = velocidad

    semilla_supersonica = 1.05 * v_c
    for i in indices_supersonicos:
        velocidad = encontrar_velocidad_por_newton_raphson(
            arreglo_r[i], v_c, r_c, semilla_supersonica,
            tolerancia=tolerancia, max_iter=max_iter
        )
        velocidades[i] = velocidad
        if np.isfinite(velocidad):
            semilla_supersonica = velocidad

    return velocidades
