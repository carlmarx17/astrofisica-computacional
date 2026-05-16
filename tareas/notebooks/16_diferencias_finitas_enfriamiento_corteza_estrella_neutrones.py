#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Tarea 16
# Enfriamiento termico de la corteza de una estrella de neutrones
# Metodo explicito FTCS con termino de perdida por neutrinos
# ============================================================

# Parametros fisicos
L = 1000.0          # profundidad de la corteza [m]
alpha = 1.0         # difusividad termica
eps_nu = 1e-45      # constante de emision de neutrinos
T_inicial = 1e9     # temperatura inicial uniforme [K]

# Condiciones de frontera de Dirichlet
T_superficie = 1e6  # temperatura en z = 0 [K]
T_interior = 1e9    # temperatura en z = L [K]

# Parametros numericos
Nz = 100
dz = L / (Nz - 1)
dt = 0.4 * dz**2 / alpha
Nt = 5000


def evolucionar_temperatura():
    """
    Resuelve la ecuacion de difusion con un termino sumidero:

        dT/dt = alpha * d2T/dz2 - eps_nu * T^5

    usando diferencias finitas explicitas (FTCS).
    """
    z = np.linspace(0.0, L, Nz)
    T = np.full(Nz, T_inicial, dtype=float)
    T[0] = T_superficie
    T[-1] = T_interior

    perfiles = {0: T.copy()}

    pasos_para_graficar = [100, 1000, 2500, 4999]

    for n in range(Nt):
        T_nueva = T.copy()

        difusion = alpha * (T[2:] - 2.0 * T[1:-1] + T[:-2]) / dz**2
        sumidero = eps_nu * T[1:-1] ** 5
        T_nueva[1:-1] = T[1:-1] + dt * (difusion - sumidero)

        T_nueva[0] = T_superficie
        T_nueva[-1] = T_interior
        T = T_nueva

        if n in pasos_para_graficar:
            perfiles[n] = T.copy()

    return z, perfiles, pasos_para_graficar


z, perfiles, pasos_para_graficar = evolucionar_temperatura()

print("=" * 64)
print("  ENFRIAMIENTO TERMICO DE LA CORTEZA DE UNA ESTRELLA DE NEUTRONES")
print("=" * 64)
print(f"  Profundidad total L          = {L:.1f} m")
print(f"  Difusividad termica alpha    = {alpha:.3e}")
print(f"  Constante de neutrinos epsnu = {eps_nu:.3e}")
print(f"  Temperatura inicial          = {T_inicial:.3e} K")
print(f"  Temperatura superficial      = {T_superficie:.3e} K")
print(f"  Temperatura interior         = {T_interior:.3e} K")
print()
print(f"  Nz = {Nz}")
print(f"  dz = {dz:.4f} m")
print(f"  dt = {dt:.4f} s")
print(f"  Nt = {Nt}")
print(f"  Tiempo total simulado = {Nt * dt:.4f} s")
print("=" * 64)


plt.figure(figsize=(10, 6))
plt.plot(z, perfiles[0], linestyle="--", linewidth=2.0, label="t = 0 s")

for paso in pasos_para_graficar:
    tiempo_actual = paso * dt
    plt.plot(z, perfiles[paso], linewidth=2.0, label=f"Paso {paso} (t ~ {tiempo_actual:.2f} s)")

plt.title("Enfriamiento termico de la corteza de una estrella de neutrones")
plt.xlabel("Profundidad z [m]")
plt.ylabel("Temperatura T [K]")
plt.grid(True, linestyle=":", alpha=0.7)
plt.legend()
plt.tight_layout()
plt.show()
