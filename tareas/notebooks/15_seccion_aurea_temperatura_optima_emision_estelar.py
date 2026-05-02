#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Parametros del modelo
# ============================================================
sigma = 5.67e-8   # W m^-2 K^-4 (constante de Stefan-Boltzmann)
T0 = 10000        # K
hv_kB = 5000      # K  (h*nu / k_B)

T_min = 3000      # K
T_max = 50000     # K


def P(T):
    """
    Potencia neta del modelo de emision estelar.
    """
    return sigma * T**4 * np.exp(-T / T0) * (1 - np.exp(-hv_kB / T))


def seccion_aurea(f, a, b, eps=50):
    """
    Maximiza f en [a, b] usando la seccion aurea.

    eps : tolerancia de convergencia en K
    """
    phi = (np.sqrt(5) - 1) / 2
    iteraciones = []

    x1 = b - phi * (b - a)
    x2 = a + phi * (b - a)
    f1 = f(x1)
    f2 = f(x2)

    while (b - a) > eps:
        iteraciones.append({
            "a": a,
            "b": b,
            "x1": x1,
            "x2": x2,
            "f(x1)": f1,
            "f(x2)": f2,
            "intervalo": b - a,
        })

        if f1 < f2:
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + phi * (b - a)
            f2 = f(x2)
        else:
            b = x2
            x2 = x1
            f2 = f1
            x1 = b - phi * (b - a)
            f1 = f(x1)

    T_opt = 0.5 * (a + b)
    P_opt = f(T_opt)

    return T_opt, P_opt, iteraciones


T_opt, P_opt, iters = seccion_aurea(P, T_min, T_max, eps=50)

print("=" * 60)
print("  TEMPERATURA OPTIMA DE EMISION ESTELAR")
print("=" * 60)
print("\n  Parametros:")
print(f"    sigma   = {sigma:.2e}  W m^-2 K^-4")
print(f"    T0      = {T0}  K")
print(f"    hnu/k_B = {hv_kB}  K")
print(f"    Intervalo: [{T_min}, {T_max}] K")
print("\n--- Seccion Aurea (eps = 50 K) ---")
print(f"  Iteraciones realizadas : {len(iters)}")
print(f"  T optima encontrada    : {T_opt:.2f} K")
print(f"  P(T_opt)               : {P_opt:.6e} W m^-2")

print(f"\n  {'Iter':>4}  {'a':>10}  {'b':>10}  {'x1':>10}  {'x2':>10}  {'Intervalo':>12}")
print("  " + "-" * 62)
for i, it in enumerate(iters[:10]):
    print(
        f"  {i + 1:>4}  {it['a']:>10.1f}  {it['b']:>10.1f}  "
        f"{it['x1']:>10.1f}  {it['x2']:>10.1f}  {it['intervalo']:>12.2f}"
    )
if len(iters) > 10:
    print(f"  ...  ({len(iters) - 10} iteraciones adicionales)")


T_vals = np.linspace(T_min, T_max, 5000)
P_vals = P(T_vals)

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(T_vals, P_vals, color="#2563EB", linewidth=2.2, label=r"$P(T)$")
ax.axvline(T_opt, color="#DC2626", linestyle="--", linewidth=1.4, alpha=0.8)
ax.scatter(
    [T_opt],
    [P_opt],
    color="#DC2626",
    s=80,
    zorder=5,
    label=f"Maximo: T = {T_opt:.0f} K\nP = {P_opt:.3e} W m^-2",
)

ax.annotate(
    f"  T* = {T_opt:.0f} K\n  P* = {P_opt:.3e} W/m^2",
    xy=(T_opt, P_opt),
    xytext=(T_opt + 5000, P_opt * 0.85),
    fontsize=10,
    arrowprops=dict(arrowstyle="->", color="#374151"),
    color="#374151",
)

ax.set_xlabel("Temperatura T [K]", fontsize=12)
ax.set_ylabel(r"Potencia neta $P(T)$ [W m$^{-2}$]", fontsize=12)
ax.set_title(
    "Temperatura Optima de Emision Estelar\n"
    r"$P(T) = \sigma T^4 e^{-T/T_0}\left(1 - e^{-h\nu/(k_B T)}\right)$",
    fontsize=13,
)
ax.legend(fontsize=11)
ax.set_xlim(T_min, T_max)
ax.yaxis.get_major_formatter().set_powerlimits((0, 0))
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("tareas/notebooks/15_grafica_P_T.png", dpi=150, bbox_inches="tight")
plt.show()

print("\n  Grafica guardada como '15_grafica_P_T.png'")
print("=" * 60)
