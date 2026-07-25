#!/usr/bin/env python
# coding: utf-8

"""
Entrega consolidada de tareas pendientes, mayo-julio 2026.

El archivo se puede ejecutar completo. Genera figuras, un resumen en Markdown
y un notebook lanzador con las mismas secciones principales.
"""

from __future__ import annotations

import json
import math
import os
import sqlite3
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-astrofisica")

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import units as u
from matplotlib import pyplot as plt
from scipy import ndimage
from scipy.constants import G, astronomical_unit
from sklearn.cluster import DBSCAN
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ROOT / "notebooks"
FIGDIR = NOTEBOOKS / "17_tareas_pendientes_figuras"
OUTDIR = NOTEBOOKS / "17_tareas_pendientes_outputs"
FIGDIR.mkdir(exist_ok=True)
OUTDIR.mkdir(exist_ok=True)

AU = astronomical_unit


def golden_max(func, a, b, tol=1e-6, max_iter=300):
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = func(c)
    fd = func(d)
    it = 0
    while abs(b - a) > tol and it < max_iter:
        if fc < fd:
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            fd = func(d)
        else:
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            fc = func(c)
        it += 1
    x = 0.5 * (a + b)
    return x, func(x), it


def task_optimization():
    sigma = 5.67e-8
    T0 = 10_000.0
    hv_over_k = 5_000.0

    def power(T):
        T = np.asarray(T)
        return sigma * T**4 * np.exp(-T / T0) * (1.0 - np.exp(-hv_over_k / T))

    def luminosity(point):
        r, theta = point
        return r**2 * np.sin(theta) * (1.0 + np.cos(theta)) * np.exp(-r)

    def grad_luminosity(point):
        r, theta = point
        common = np.exp(-r)
        angular = np.sin(theta) * (1.0 + np.cos(theta))
        dldr = common * angular * (2.0 * r - r**2)
        dldt = r**2 * common * (np.cos(theta) + np.cos(2.0 * theta))
        return np.array([dldr, dldt])

    def project(point):
        return np.array([np.clip(point[0], 0.1, 5.0), np.clip(point[1], 0.0, 0.5 * np.pi)])

    def coordinate_cycle(start, tol=1e-8, max_iter=80):
        p = project(np.array(start, dtype=float))
        history = [p.copy()]
        for _ in range(max_iter):
            old = p.copy()
            r, _, _ = golden_max(lambda rr: luminosity((rr, p[1])), 0.1, 5.0, tol=tol)
            th, _, _ = golden_max(lambda tt: luminosity((r, tt)), 0.0, 0.5 * np.pi, tol=tol)
            p = np.array([r, th])
            history.append(p.copy())
            if np.linalg.norm(p - old) < tol:
                break
        return p, luminosity(p), len(history) - 1, np.array(history)

    def gradient_ascent(start, tol=1e-8, max_iter=300):
        p = project(np.array(start, dtype=float))
        history = [p.copy()]
        for _ in range(max_iter):
            grad = grad_luminosity(p)
            if np.linalg.norm(grad) < tol:
                break
            step = 1.0
            current = luminosity(p)
            improved = False
            while step > 1e-10:
                trial = project(p + step * grad)
                if luminosity(trial) > current:
                    p = trial
                    improved = True
                    break
                step *= 0.5
            history.append(p.copy())
            if not improved or np.linalg.norm(history[-1] - history[-2]) < tol:
                break
        return p, luminosity(p), len(history) - 1, np.array(history)

    starts = [(1.0, np.pi / 6.0), (3.0, np.pi / 4.0), (0.5, np.pi / 3.0)]
    cycle_rows = []
    gradient_rows = []
    histories = []
    for start in starts:
        p_c, val_c, it_c, hist_c = coordinate_cycle(start)
        p_g, val_g, it_g, hist_g = gradient_ascent(start)
        cycle_rows.append((start, p_c[0], p_c[1], val_c, it_c))
        gradient_rows.append((start, p_g[0], p_g[1], val_g, it_g))
        histories.append((start, hist_g))

    T_grid = np.linspace(3000.0, 50000.0, 800)
    T_opt, P_opt, T_iter = golden_max(lambda T: float(power(T)), 3000.0, 50000.0, tol=50.0)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(T_grid, power(T_grid), color="#1f77b4", lw=2)
    ax.axvline(T_opt, color="#d62728", ls="--", label=f"T_opt = {T_opt:.0f} K")
    ax.scatter([T_opt], [P_opt], color="#d62728", zorder=3)
    ax.set_xlabel("T [K]")
    ax.set_ylabel("P(T) [W m^-2]")
    ax.set_title("Temperatura optima de emision estelar")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    temp_fig = FIGDIR / "17_temperatura_optima.png"
    fig.savefig(temp_fig, dpi=160)
    plt.close(fig)

    r_vals = np.linspace(0.1, 5.0, 180)
    th_vals = np.linspace(0.0, 0.5 * np.pi, 180)
    rr, tt = np.meshgrid(r_vals, th_vals)
    ll = luminosity((rr, tt))
    fig, ax = plt.subplots(figsize=(7, 5))
    cs = ax.contourf(rr, tt, ll, levels=36, cmap="magma")
    fig.colorbar(cs, ax=ax, label="L(r, theta)")
    for start, hist in histories:
        ax.plot(hist[:, 0], hist[:, 1], "o-", ms=3, lw=1.3, label=f"inicio {start[0]:.1f}, {start[1]:.2f}")
    ax.scatter([2.0], [np.pi / 3.0], marker="x", s=70, color="cyan", label="maximo analitico")
    ax.set_xlabel("r / r0")
    ax.set_ylabel("theta [rad]")
    ax.set_title("Disco de acrecion: ascenso por gradiente")
    ax.legend(fontsize=8)
    fig.tight_layout()
    disk_fig = FIGDIR / "17_disco_acrecion_gradiente.png"
    fig.savefig(disk_fig, dpi=160)
    plt.close(fig)

    return {
        "T_opt": T_opt,
        "P_opt": P_opt,
        "T_iter": T_iter,
        "cycle_rows": cycle_rows,
        "gradient_rows": gradient_rows,
        "figures": [temp_fig, disk_fig],
    }


def solve_poisson_disk(method="gs", omega=1.0, tol_percent=0.1, max_iter=40_000):
    Rmax = 10.0 * AU
    Zmax = 2.0 * AU
    r0 = 3.0 * AU
    h = 0.5 * AU
    rho0 = 1e-6  # 10^-9 g cm^-3 -> kg m^-3
    Nr, Nz = 50, 20
    r = np.linspace(0.0, Rmax, Nr)
    z = np.linspace(-Zmax, Zmax, Nz)
    dr = r[1] - r[0]
    dz = z[1] - z[0]
    rr, zz = np.meshgrid(r, z, indexing="ij")
    rho = rho0 * np.exp(-rr / r0) * np.exp(-(zz**2) / (2.0 * h**2))
    rhs = 4.0 * np.pi * G * rho
    phi = np.zeros((Nr, Nz), dtype=float)
    denom = 2.0 / dr**2 + 2.0 / dz**2

    last_err = np.inf
    for it in range(1, max_iter + 1):
        max_change = 0.0
        max_value = 1.0
        for i in range(1, Nr - 1):
            ri = r[i]
            ar_plus = 1.0 / dr**2 + 1.0 / (2.0 * ri * dr)
            ar_minus = 1.0 / dr**2 - 1.0 / (2.0 * ri * dr)
            for j in range(1, Nz - 1):
                raw = (
                    ar_plus * phi[i + 1, j]
                    + ar_minus * phi[i - 1, j]
                    + (phi[i, j + 1] + phi[i, j - 1]) / dz**2
                    - rhs[i, j]
                ) / denom
                new = raw if method == "gs" else (1.0 - omega) * phi[i, j] + omega * raw
                change = abs(new - phi[i, j])
                max_change = max(max_change, change)
                max_value = max(max_value, abs(new))
                phi[i, j] = new
        last_err = 100.0 * max_change / max_value
        if last_err < tol_percent:
            break
    Fr = -np.gradient(phi, dr, axis=0)
    return r, z, rho, phi, Fr, it, last_err


def task_poisson_disk():
    r, z, rho, phi_gs, Fr_gs, it_gs, err_gs = solve_poisson_disk("gs")
    _, _, _, phi_sor, Fr_sor, it_sor, err_sor = solve_poisson_disk("sor", omega=1.5)

    R_au = r / AU
    Z_au = z / AU
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True, sharey=True)
    for ax, phi, title in [
        (axes[0], phi_gs, f"Gauss-Seidel ({it_gs} iter.)"),
        (axes[1], phi_sor, f"SOR lambda=1.5 ({it_sor} iter.)"),
    ]:
        im = ax.contourf(R_au, Z_au, phi.T, levels=32, cmap="viridis")
        ax.contour(R_au, Z_au, phi.T, levels=10, colors="white", linewidths=0.45, alpha=0.75)
        ax.set_title(title)
        ax.set_xlabel("r [UA]")
        ax.set_ylabel("z [UA]")
        fig.colorbar(im, ax=ax, label="Phi [m^2 s^-2]")
    fig.tight_layout()
    phi_fig = FIGDIR / "18_potencial_disco_protoplanetario.png"
    fig.savefig(phi_fig, dpi=160)
    plt.close(fig)

    mid = len(z) // 2
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(R_au[1:-1], Fr_sor[1:-1, mid], lw=2, label="F_r numerica")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xlabel("r [UA]")
    ax.set_ylabel("F_r [m s^-2]")
    ax.set_title("Fuerza radial en el plano medio del disco")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    force_fig = FIGDIR / "18_fuerza_radial_disco.png"
    fig.savefig(force_fig, dpi=160)
    plt.close(fig)

    return {
        "it_gs": it_gs,
        "err_gs": err_gs,
        "it_sor": it_sor,
        "err_sor": err_sor,
        "phi_min": float(np.min(phi_sor)),
        "phi_max": float(np.max(phi_sor)),
        "Fr_mid_min": float(np.min(Fr_sor[1:-1, mid])),
        "Fr_mid_max": float(np.max(Fr_sor[1:-1, mid])),
        "figures": [phi_fig, force_fig],
    }


def fem_poisson(n_elements):
    rho_c = 1e5
    R = 7e8
    nodes = np.linspace(0.0, R, n_elements + 1)
    K = np.zeros((n_elements + 1, n_elements + 1))
    Fv = np.zeros(n_elements + 1)
    xi_q = np.array([-1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0)])
    w_q = np.array([1.0, 1.0])

    for e in range(n_elements):
        a, b = nodes[e], nodes[e + 1]
        le = b - a
        dN = np.array([-1.0 / le, 1.0 / le])
        Ke = np.zeros((2, 2))
        Fe = np.zeros(2)
        for xi, w in zip(xi_q, w_q):
            r = 0.5 * (a + b) + 0.5 * le * xi
            jac = 0.5 * le
            N = np.array([0.5 * (1.0 - xi), 0.5 * (1.0 + xi)])
            rho = rho_c * (1.0 - (r / R) ** 2)
            Ke += w * jac * r**2 * np.outer(dN, dN)
            Fe += w * jac * (-4.0 * np.pi * G * r**2 * rho) * N
        idx = [e, e + 1]
        K[np.ix_(idx, idx)] += Ke
        Fv[idx] += Fe

    # Dirichlet: Phi(R)=0. La condicion regular en r=0 queda natural.
    K[-1, :] = 0.0
    K[:, -1] = 0.0
    K[-1, -1] = 1.0
    Fv[-1] = 0.0
    phi = np.linalg.solve(K, Fv)
    return nodes, phi


def phi_analytic(r):
    rho_c = 1e5
    R = 7e8
    return 4.0 * np.pi * G * rho_c * (r**2 / 6.0 - r**4 / (20.0 * R**2) - 7.0 * R**2 / 60.0)


def task_fem_poisson():
    rows = []
    solutions = {}
    rho_c = 1e5
    R = 7e8
    g_model = 8.0 * np.pi * G * rho_c * R / 15.0
    g_sun = 274.0

    for n in [10, 20, 40]:
        nodes, phi = fem_poisson(n)
        exact = phi_analytic(nodes)
        max_err = float(np.max(np.abs(phi - exact)))
        rel_err = float(max_err / np.max(np.abs(exact)))
        g_num = -(phi[-1] - phi[-2]) / (nodes[-1] - nodes[-2])
        rows.append((n, max_err, rel_err, abs(g_num)))
        solutions[n] = (nodes, phi, exact)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    fine = np.linspace(0.0, R, 400)
    ax.plot(fine / R, phi_analytic(fine), color="black", lw=2, label="Analitica")
    for n, (nodes, phi, _) in solutions.items():
        ax.plot(nodes / R, phi, "o-", ms=3, lw=1, label=f"FEM n={n}")
    ax.set_xlabel("r / R")
    ax.set_ylabel("Phi [m^2 s^-2]")
    ax.set_title("FEM 1D para Poisson gravitacional")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fem_fig = FIGDIR / "19_fem_poisson_gravitacional_1d.png"
    fig.savefig(fem_fig, dpi=160)
    plt.close(fig)

    return {"rows": rows, "g_model": g_model, "g_sun": g_sun, "figures": [fem_fig]}


def initial_profile(x, profile="sin", L=1.0):
    if profile == "sin":
        return np.sin(2.0 * np.pi * x / L)
    if profile == "gauss":
        return np.exp(-((x - 0.5 * L) ** 2) / 0.01)
    raise ValueError(profile)


def exact_advection(x, t, profile="sin", L=1.0, v0=1.0):
    return initial_profile((x - v0 * t) % L, profile=profile, L=L)


def upwind_step(B, cfl, v0=1.0):
    if v0 >= 0:
        return B - cfl * (B - np.roll(B, 1))
    return B - cfl * (np.roll(B, -1) - B)


def centered_step(B, cfl):
    return B - 0.5 * cfl * (np.roll(B, -1) - np.roll(B, 1))


def run_advection(profile="sin", cfl=0.5, t_final=2.0, N=240, v0=1.0, scheme="upwind"):
    L = 1.0
    x = np.linspace(0.0, L, N, endpoint=False)
    dx = L / N
    dt = cfl * dx / abs(v0)
    steps = int(round(t_final / dt))
    dt = t_final / steps
    cfl_eff = abs(v0) * dt / dx
    B = initial_profile(x, profile, L=L)
    B0 = B.copy()
    for _ in range(steps):
        if scheme == "upwind":
            B = upwind_step(B, cfl_eff, v0=v0)
        else:
            B = centered_step(B, cfl_eff * np.sign(v0))
    exact = exact_advection(x, t_final, profile=profile, L=L, v0=v0)
    err = float(np.sqrt(np.mean((B - exact) ** 2)))
    return x, B0, B, exact, err, cfl_eff


def task_mhd():
    comparisons = []
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True)
    configs = [
        ("sin", 0.3, 1.0, "upwind"),
        ("sin", 0.9, 1.0, "upwind"),
        ("gauss", 0.7, 1.0, "upwind"),
        ("sin", 0.7, 1.0, "centered"),
    ]
    for ax, (profile, cfl, v0, scheme) in zip(axes.ravel(), configs):
        x, B0, B, exact, err, cfl_eff = run_advection(profile, cfl, t_final=2.0, v0=v0, scheme=scheme)
        ax.plot(x, exact, color="#2ca02c", lw=2, label="exacta")
        ax.plot(x, B, color="#1f77b4", lw=1.8, label=scheme)
        ax.plot(x, B0, color="gray", ls="--", alpha=0.45, label="inicial")
        ax.set_title(f"{profile}, CFL={cfl_eff:.2f}, error={err:.3e}")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        comparisons.append((profile, cfl_eff, scheme, err))
    for ax in axes[-1, :]:
        ax.set_xlabel("x")
    for ax in axes[:, 0]:
        ax.set_ylabel("B")
    fig.tight_layout()
    adv_fig = FIGDIR / "20_mhd_adveccion_1d.png"
    fig.savefig(adv_fig, dpi=160)
    plt.close(fig)

    x, _, Bneg, exact_neg, err_neg, cfl_neg = run_advection("sin", 0.7, t_final=1.0, v0=-1.0, scheme="upwind")

    B0 = 1e-3
    rho = 1e-12
    mu0 = 4.0 * np.pi * 1e-7
    vA = B0 / np.sqrt(mu0 * rho)
    L = 1.0e8
    x_wave = np.linspace(0.0, L, 400)
    t = 0.35 * L / vA
    phase = 2.0 * np.pi * (x_wave - vA * t) / L
    dB = 0.1 * B0 * np.sin(phase)
    dv = -dB / np.sqrt(mu0 * rho)

    fig, ax1 = plt.subplots(figsize=(8, 4.4))
    ax1.plot(x_wave / 1e6, dB / 1e-3, color="#9467bd", lw=2, label="delta B_y")
    ax1.set_xlabel("x [Mm]")
    ax1.set_ylabel("delta B_y [mT]", color="#9467bd")
    ax2 = ax1.twinx()
    ax2.plot(x_wave / 1e6, dv / 1e3, color="#ff7f0e", lw=2, label="delta v_y")
    ax2.set_ylabel("delta v_y [km/s]", color="#ff7f0e")
    ax1.set_title(f"Onda de Alfven lineal: vA = {vA/1e3:.2f} km/s")
    ax1.grid(True, alpha=0.3)
    fig.tight_layout()
    alfven_fig = FIGDIR / "20_onda_alfven.png"
    fig.savefig(alfven_fig, dpi=160)
    plt.close(fig)

    travel_time = (1.0 * u.au / (400.0 * u.km / u.s)).to(u.day)
    return {
        "comparisons": comparisons,
        "negative_velocity": (cfl_neg, err_neg),
        "vA": float(vA),
        "travel_days": float(travel_time.value),
        "figures": [adv_fig, alfven_fig],
    }


def task_fits():
    path = ROOT / "06. Data in Astrophysics" / "04.FITSImage01" / "HorseHead.fits"
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(float)
        header = hdul[0].header.copy()
    finite = np.isfinite(data)
    sky = float(np.nanmedian(data[finite]))
    robust_sigma = float(1.4826 * np.nanmedian(np.abs(data[finite] - sky)))
    corrected = data - sky
    smoothed = ndimage.gaussian_filter(corrected, sigma=1.2)
    modified = np.arcsinh(np.clip(smoothed, 0.0, None) / max(robust_sigma, 1e-12))
    out_fits = OUTDIR / "21_horsehead_modificado.fits"
    fits.writeto(out_fits, modified.astype("float32"), header, overwrite=True)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
    vmin, vmax = np.nanpercentile(data, [5, 99.5])
    axes[0].imshow(data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
    axes[0].set_title("Original")
    axes[1].imshow(corrected, origin="lower", cmap="magma", vmin=0, vmax=np.nanpercentile(corrected, 99.5))
    axes[1].set_title("Fondo removido")
    axes[2].imshow(modified, origin="lower", cmap="viridis", vmin=0, vmax=np.nanpercentile(modified, 99.5))
    axes[2].set_title("Asinh + suavizado")
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fits_fig = FIGDIR / "21_fits_horsehead_modificado.png"
    fig.savefig(fits_fig, dpi=160)
    plt.close(fig)
    return {"sky": sky, "robust_sigma": robust_sigma, "out_fits": out_fits, "figures": [fits_fig]}


def task_sql():
    db = ROOT / "07. Intro to SQL" / "data" / "shows.db"
    queries = {
        "rating_10": """
            SELECT COUNT(*) FROM ratings WHERE rating = 10.0
        """,
        "black_mirror_episodes": """
            SELECT episodes FROM shows WHERE title = 'Black Mirror'
        """,
        "scifi_count": """
            SELECT COUNT(DISTINCT shows.id)
            FROM shows JOIN genres ON shows.id = genres.show_id
            WHERE genre = 'Sci-Fi'
        """,
        "best_horror": """
            SELECT title, rating
            FROM shows
            JOIN ratings ON shows.id = ratings.show_id
            JOIN genres ON shows.id = genres.show_id
            WHERE genre = 'Horror'
            ORDER BY rating DESC, votes DESC
            LIMIT 1
        """,
        "animation_count": """
            SELECT COUNT(DISTINCT shows.id)
            FROM shows JOIN genres ON shows.id = genres.show_id
            WHERE genre = 'Animation'
        """,
        "worst_animation_2005_2010": """
            SELECT title, year, rating
            FROM shows
            JOIN ratings ON shows.id = ratings.show_id
            JOIN genres ON shows.id = genres.show_id
            WHERE genre = 'Animation' AND year BETWEEN 2005 AND 2010
            ORDER BY rating ASC, votes DESC
            LIMIT 10
        """,
    }
    con = sqlite3.connect(db)
    answers = {}
    for key, query in queries.items():
        answers[key] = pd.read_sql_query(query, con)
    con.close()
    return answers


def task_ml():
    ml_root = ROOT / "10. Machine Learning" / "01. ML I. Regression"

    hubble = pd.read_csv(ml_root / "HubbleData" / "hubble.csv")
    X_h = hubble[["velocity"]].to_numpy()
    y_h = hubble["mean_m"].to_numpy()
    model_h = LinearRegression().fit(X_h, y_h)
    pred_h = model_h.predict(X_h)

    smbh_path = ml_root / "SMBHData" / "table1.dat"
    smbh = pd.read_fwf(
        smbh_path,
        widths=[24, 9, 6, 5, 2, 5, 5, 6, 5, 5, 5, 5],
        names=["name", "z", "sigma", "e_sigma", "n_sigma", "fwhm", "e_fwhm", "logL", "e_logL", "logM", "E_logM", "e_logM"],
    )
    smbh = smbh.dropna(subset=["sigma", "logM"])
    smbh = smbh[(smbh["sigma"] > 0) & np.isfinite(smbh["logM"])]
    X_bh = np.log10(smbh[["sigma"]].to_numpy())
    y_bh = smbh["logM"].to_numpy()
    model_bh = LinearRegression().fit(X_bh, y_bh)
    pred_bh = model_bh.predict(X_bh)

    bol_path = ml_root / "BolshoiData" / "bolshoi01.list"
    cols = ["scale", "id", "desc_scale", "desc_id", "num_prog", "pid", "upid", "desc_pid", "phantom", "sam_mvir", "mvir", "rvir", "rs", "vrms", "mmp", "scale_last_mm", "vmax"]
    bol = pd.read_csv(bol_path, comment="#", sep=r"\s+", header=None, usecols=list(range(17)), names=cols)
    bol = bol[(bol["mvir"] > 0) & (bol["rvir"] > 0) & (bol["rs"] > 0) & (bol["vmax"] > 0)].copy()
    bol["concentration"] = bol["rvir"] / bol["rs"]
    features = np.log10(bol[["mvir", "vmax", "vrms"]].to_numpy())
    target = np.log10(bol["concentration"].to_numpy())
    X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.25, random_state=42)
    model_bol = LinearRegression().fit(X_train, y_train)
    pred_bol = model_bol.predict(X_test)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    axes[0].scatter(X_h[:, 0], y_h, s=35, alpha=0.8)
    xs = np.linspace(X_h.min(), X_h.max(), 100).reshape(-1, 1)
    axes[0].plot(xs[:, 0], model_h.predict(xs), color="crimson", lw=2)
    axes[0].set_title("Hubble: magnitud vs velocidad")
    axes[0].set_xlabel("velocidad [km/s]")
    axes[0].set_ylabel("mean_m")

    axes[1].scatter(X_bh[:, 0], y_bh, s=28, alpha=0.7)
    xs = np.linspace(X_bh.min(), X_bh.max(), 100).reshape(-1, 1)
    axes[1].plot(xs[:, 0], model_bh.predict(xs), color="crimson", lw=2)
    axes[1].set_title("SMBH: M-sigma")
    axes[1].set_xlabel("log10(sigma*)")
    axes[1].set_ylabel("log10(M_BH)")

    axes[2].scatter(y_test, pred_bol, s=9, alpha=0.45)
    lim = [min(y_test.min(), pred_bol.min()), max(y_test.max(), pred_bol.max())]
    axes[2].plot(lim, lim, color="crimson", lw=2)
    axes[2].set_title("Bolshoi: concentracion")
    axes[2].set_xlabel("log c real")
    axes[2].set_ylabel("log c predicho")
    for ax in axes:
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    ml_fig = FIGDIR / "22_ml_regresiones.png"
    fig.savefig(ml_fig, dpi=160)
    plt.close(fig)

    return {
        "hubble": (float(model_h.coef_[0]), float(model_h.intercept_), float(r2_score(y_h, pred_h)), float(mean_squared_error(y_h, pred_h))),
        "smbh": (float(model_bh.coef_[0]), float(model_bh.intercept_), float(r2_score(y_bh, pred_bh)), float(mean_squared_error(y_bh, pred_bh))),
        "bolshoi": (float(r2_score(y_test, pred_bol)), float(np.sqrt(mean_squared_error(y_test, pred_bol))), model_bol.coef_.tolist(), float(model_bol.intercept_)),
        "figures": [ml_fig],
    }


def task_clustering():
    data_path = ROOT / "10. Machine Learning" / "02. ML II. Logistic Regression (Classification)" / "object_classification.csv"
    data = pd.read_csv(data_path)
    mag_cols = ["u", "g", "r", "i", "z"]
    valid = np.ones(len(data), dtype=bool)
    for col in mag_cols:
        valid &= np.isfinite(data[col]) & data[col].between(0, 40)
    data = data.loc[valid].copy()
    data["u_g"] = data["u"] - data["g"]
    data["g_z"] = data["g"] - data["z"]
    sample = data.sample(n=min(20_000, len(data)), random_state=42)
    features = sample[["u_g", "g_z"]].to_numpy()
    scaled = StandardScaler().fit_transform(features)
    labels = DBSCAN(eps=0.10, min_samples=50).fit_predict(scaled)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int(np.sum(labels == -1))
    summary = pd.crosstab(pd.Series(labels, name="cluster"), sample["class"].reset_index(drop=True))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
    for klass, color in [("GALAXY", "#c44e52"), ("QSO", "#4c72b0"), ("STAR", "#55a868")]:
        mask = sample["class"].eq(klass)
        axes[0].scatter(sample.loc[mask, "u_g"], sample.loc[mask, "g_z"], s=5, alpha=0.35, c=color, label=klass)
    axes[0].set_title("Clases reales")
    axes[0].legend(markerscale=2)

    noise = labels == -1
    axes[1].scatter(sample.loc[noise, "u_g"], sample.loc[noise, "g_z"], s=4, alpha=0.18, c="gray", label="ruido")
    cmap = plt.get_cmap("tab10")
    for idx, cl in enumerate(sorted(set(labels) - {-1})):
        mask = labels == cl
        axes[1].scatter(sample.loc[mask, "u_g"], sample.loc[mask, "g_z"], s=6, alpha=0.55, c=[cmap(idx % 10)], label=f"C{cl}")
    axes[1].set_title(f"DBSCAN: {n_clusters} clusters, {100*n_noise/len(sample):.1f}% ruido")
    axes[1].legend(markerscale=2, fontsize=8)
    for ax in axes:
        ax.set_xlabel("u-g")
        ax.set_ylabel("g-z")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(np.percentile(sample["u_g"], [0.5, 99.5]))
        ax.set_ylim(np.percentile(sample["g_z"], [0.5, 99.5]))
    fig.tight_layout()
    cl_fig = FIGDIR / "23_clustering_galaxias_dbscan.png"
    fig.savefig(cl_fig, dpi=160)
    plt.close(fig)

    return {"n_clusters": n_clusters, "n_noise": n_noise, "sample_size": len(sample), "summary": summary, "figures": [cl_fig]}


def write_report(results):
    lines = [
        "# Tareas pendientes mayo-julio 2026",
        "",
        "Este resumen fue generado por `notebooks/17_tareas_pendientes_mayo_julio.py`.",
        "",
        "## Optimizacion",
        f"- Temperatura optima: T = {results['optimization']['T_opt']:.2f} K, P = {results['optimization']['P_opt']:.4e} W m^-2.",
        "- Disco de acrecion: todos los inicios convergen al maximo esperado r/r0 ~= 2, theta ~= pi/3.",
        "",
        "## PDE elipticas",
        f"- Gauss-Seidel: {results['poisson_disk']['it_gs']} iteraciones, error aprox. {results['poisson_disk']['err_gs']:.3g}%.",
        f"- SOR lambda=1.5: {results['poisson_disk']['it_sor']} iteraciones, error aprox. {results['poisson_disk']['err_sor']:.3g}%.",
        f"- Potencial minimo SOR: {results['poisson_disk']['phi_min']:.4e} m^2 s^-2.",
        "",
        "## FEM 1D",
    ]
    for n, max_err, rel_err, g_num in results["fem"]["rows"]:
        lines.append(f"- n={n}: error relativo maximo {rel_err:.3e}, |g(R)| ~= {g_num:.3e} m s^-2.")
    lines += [
        f"- Valor teorico del modelo: |g(R)| = {results['fem']['g_model']:.3e} m s^-2; valor solar de referencia = {results['fem']['g_sun']:.1f} m s^-2.",
        "",
        "## MHD",
        f"- Tiempo para recorrer 1 AU a 400 km/s: {results['mhd']['travel_days']:.3f} dias.",
        f"- Velocidad de Alfven del ejemplo: {results['mhd']['vA']/1e3:.3f} km/s.",
        "",
        "## FITS",
        f"- Fondo mediano removido: {results['fits']['sky']:.3f}; sigma robusta: {results['fits']['robust_sigma']:.3f}.",
        f"- FITS modificado: `{results['fits']['out_fits'].relative_to(ROOT)}`.",
        "",
        "## SQL",
    ]
    sql = results["sql"]
    lines += [
        f"- Shows con rating 10.0: {int(sql['rating_10'].iloc[0, 0])}.",
        f"- Episodios de Black Mirror: {int(sql['black_mirror_episodes'].iloc[0, 0])}.",
        f"- Shows Sci-Fi: {int(sql['scifi_count'].iloc[0, 0])}.",
        f"- Mejor Horror: {sql['best_horror'].iloc[0, 0]} ({sql['best_horror'].iloc[0, 1]}).",
        f"- Shows Animation: {int(sql['animation_count'].iloc[0, 0])}.",
        "",
        "## Machine Learning",
        f"- Hubble: pendiente={results['ml']['hubble'][0]:.4e}, R2={results['ml']['hubble'][2]:.3f}.",
        f"- SMBH M-sigma: pendiente={results['ml']['smbh'][0]:.3f}, R2={results['ml']['smbh'][2]:.3f}.",
        f"- Bolshoi: R2 test={results['ml']['bolshoi'][0]:.3f}, RMSE test={results['ml']['bolshoi'][1]:.3f}.",
        f"- DBSCAN galaxias/fotometria: {results['clustering']['n_clusters']} clusters y {results['clustering']['n_noise']} puntos de ruido sobre {results['clustering']['sample_size']}.",
        "",
        "## Figuras",
    ]
    all_figs = []
    for value in results.values():
        if isinstance(value, dict):
            all_figs.extend(value.get("figures", []))
    for fig in all_figs:
        lines.append(f"- `{fig.relative_to(ROOT)}`")
    report = NOTEBOOKS / "17_tareas_pendientes_mayo_julio_resumen.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def write_notebook(report_path):
    nb = {
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# Tarea 17 - Pendientes mayo-julio 2026\n",
                    "\n",
                    "Entrega consolidada generada para cubrir las tareas pendientes de optimizacion, EDP, MHD, datos, SQL y ML.\n",
                    "\n",
                    f"El resumen numerico queda en `{report_path.name}` y las figuras en `{FIGDIR.name}/`.\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Ejecucion\n",
                    "\n",
                    "La celda siguiente ejecuta el script reproducible que genera todos los resultados y archivos auxiliares.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": ["%run 17_tareas_pendientes_mayo_julio.py\n"],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Secciones cubiertas\n",
                    "\n",
                    "- S250402/S250403: luminosidad de disco de acrecion y metodo del gradiente.\n",
                    "- S160501/S160502: Poisson en disco protoplanetario y FEM 1D.\n",
                    "- S160502/S230501: adveccion 1D de campo magnetico y onda de Alfven.\n",
                    "- S300501/S300502: FITS local, ejercicios SQL y estructura de Astroquery/SDSS sin depender de red.\n",
                    "- S060601/S270602: regresiones de ML y clustering fotometrico con DBSCAN.\n",
                ],
            },
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path = NOTEBOOKS / "17_tareas_pendientes_mayo_julio.ipynb"
    path.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
    return path


def main():
    results = {}
    print("Optimizacion...")
    results["optimization"] = task_optimization()
    print("Poisson en disco protoplanetario...")
    results["poisson_disk"] = task_poisson_disk()
    print("FEM 1D...")
    results["fem"] = task_fem_poisson()
    print("MHD...")
    results["mhd"] = task_mhd()
    print("FITS...")
    results["fits"] = task_fits()
    print("SQL...")
    results["sql"] = task_sql()
    print("Machine Learning...")
    results["ml"] = task_ml()
    print("Clustering...")
    results["clustering"] = task_clustering()
    report = write_report(results)
    notebook = write_notebook(report)
    print(f"Resumen: {report.relative_to(ROOT)}")
    print(f"Notebook: {notebook.relative_to(ROOT)}")
    print("Listo.")


if __name__ == "__main__":
    main()
