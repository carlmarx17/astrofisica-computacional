# =============================================================================
# constantes.py  —  El vecindario de los números
# =============================================================================
# Aquí viven todos los números físicos que necesitamos.
# Están en unidades SI (metros, kilogramos, segundos, Kelvin) para todo.
# Si cambias algo aquí, el resto del proyecto lo sigue automáticamente.
# =============================================================================

import numpy as np

# ── Constantes universales ────────────────────────────────────────────────────

GRAVEDAD          = 6.67430e-11    # Constante gravitacional de Newton  [m³ kg⁻¹ s⁻²]
BOLTZMANN         = 1.380649e-23   # Constante de Boltzmann             [J K⁻¹]
MASA_PROTON       = 1.6726219e-27  # Masa de un protón                  [kg]

# ── Parámetros solares ────────────────────────────────────────────────────────

MASA_SOL          = 1.989e30       # Masa total del Sol                 [kg]
RADIO_SOL         = 6.96e8         # Radio del Sol                      [m]
PESO_MOLECULAR    = 0.5            # µ para hidrógeno totalmente ionizado
                                   # (mitad protones, mitad electrones → µ = 0.5)

# ── Parámetros por defecto de la simulación ───────────────────────────────────

TEMPERATURA_BASE  = 1.0e6          # Temperatura típica de la corona, ¡ardiente! [K]

# Tasa de pérdida de masa del Sol. ~2×10⁻¹⁴ masas solares por año.
# Convertido: 2e-14 * MASA_SOL / (365.25 * 24 * 3600) ≈ 1.26e9 kg/s
FLUJO_MASA_BASE   = 1.26e9         # Tasa de pérdida de masa solar      [kg/s]


# ── Funciones de apoyo ────────────────────────────────────────────────────────

def convertir_flujo_masa_a_si(flujo_en_masas_solares_por_anio):
    """
    Convierte una tasa de pérdida de masa desde la unidad astronómica
    "masas solares por año" a kilogramos por segundo (SI).

    Ejemplo
    -------
    >>> convertir_flujo_masa_a_si(2e-14)   # → ~1.26e9 kg/s
    """
    segundos_por_anio = 365.25 * 24 * 3600
    return flujo_en_masas_solares_por_anio * MASA_SOL / segundos_por_anio
