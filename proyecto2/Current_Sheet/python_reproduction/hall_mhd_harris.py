"""
hall_mhd_harris.py

Standalone Python Hall-MHD evolution for the 2D Harris current sheet.

The PLUTO run in this project remains the high-resolution reference solution.
This script is a self-contained finite-difference solver that evolves the same
initial condition in Python, writes snapshots, computes diagnostics and creates
figures without reading PLUTO output.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path

MPLCONFIG = Path("/tmp") / "matplotlib-cache"
MPLCONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import SymLogNorm

# np.trapezoid exists from NumPy 2.0 onward; NumPy 1.x keeps np.trapz.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "output"

CLIMS = {
    r"$v_x$": (-0.12, 0.12),
    r"$v_y$": (-0.12, 0.12),
    r"$B_x$": (-1.0, 1.0),
    r"$B_y$": (-0.08, 0.08),
    r"$B_z$": (-0.08, 0.08),
    r"$J_z$": (-2.0, 2.0),
    r"$\nabla\cdot B$": (-1.0e-12, 1.0e-12),
}


@dataclass(frozen=True)
class Params:
    lx: float = 25.6
    ly: float = 12.8
    nx: int = 64
    ny: int = 32
    tstop: float = 5.0
    output_dt: float = 1.0
    cfl: float = 0.22
    cfl_hall: float = 0.18
    cs2: float = 0.5
    b0: float = 1.0
    width: float = 0.5
    psi0: float = 0.02
    hall_coeff: float = 1.0
    eta: float = 5.0e-3
    nu: float = 1.0e-3
    eta_h: float = 1.0e-4
    nu_h: float = 5.0e-5
    rho_floor: float = 0.05
    max_steps: int = 1_000_000
    output_dir: Path = OUT_DIR


@dataclass
class State:
    rho: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    az: np.ndarray
    bz: np.ndarray


def make_grid(p: Params) -> tuple[np.ndarray, np.ndarray, float, float]:
    dx = p.lx / p.nx
    dy = p.ly / p.ny
    x = np.linspace(-0.5 * p.lx + 0.5 * dx, 0.5 * p.lx - 0.5 * dx, p.nx)
    y = np.linspace(-0.5 * p.ly + 0.5 * dy, 0.5 * p.ly - 0.5 * dy, p.ny)
    return x, y, dx, dy


def ddx(a: np.ndarray, dx: float) -> np.ndarray:
    return (np.roll(a, -1, axis=0) - np.roll(a, 1, axis=0)) / (2.0 * dx)


def d2dx(a: np.ndarray, dx: float) -> np.ndarray:
    return (np.roll(a, -1, axis=0) - 2.0 * a + np.roll(a, 1, axis=0)) / (dx * dx)


def _pad_y_linear(a: np.ndarray, n: int) -> np.ndarray:
    """Pad y with linear ghost columns instead of flat edge repetition.

    Most variables are flat at the reflective walls after apply_boundaries().
    Az is different: its wall slope represents the saturated Bx field. Linear
    extrapolation keeps that slope and avoids a spurious current layer at the
    top and bottom boundaries.
    """
    left = a[:, :1]
    right = a[:, -1:]
    slope_left = left - a[:, 1:2]
    slope_right = right - a[:, -2:-1]
    ghosts_left = [left + k * slope_left for k in range(n, 0, -1)]
    ghosts_right = [right + k * slope_right for k in range(1, n + 1)]
    return np.concatenate(ghosts_left + [a] + ghosts_right, axis=1)


def ddy(a: np.ndarray, dy: float) -> np.ndarray:
    g = _pad_y_linear(a, 1)
    return (g[:, 2:] - g[:, :-2]) / (2.0 * dy)


def d2dy(a: np.ndarray, dy: float) -> np.ndarray:
    g = _pad_y_linear(a, 1)
    return (g[:, 2:] - 2.0 * g[:, 1:-1] + g[:, :-2]) / (dy * dy)


def lap(a: np.ndarray, dx: float, dy: float) -> np.ndarray:
    return d2dx(a, dx) + d2dy(a, dy)


def d4dx(a: np.ndarray, dx: float) -> np.ndarray:
    return (np.roll(a, -2, axis=0) - 4.0 * np.roll(a, -1, axis=0) + 6.0 * a
            - 4.0 * np.roll(a, 1, axis=0) + np.roll(a, 2, axis=0)) / dx**4


def d4dy(a: np.ndarray, dy: float) -> np.ndarray:
    g = _pad_y_linear(a, 2)
    return (g[:, 4:] - 4.0 * g[:, 3:-1] + 6.0 * g[:, 2:-2]
            - 4.0 * g[:, 1:-3] + g[:, :-4]) / dy**4


def hyp(a: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """4th-order hyper-dissipation operator: -(d⁴/dx⁴ + d⁴/dy⁴)"""
    return d4dx(a, dx) + d4dy(a, dy)


def apply_boundaries(s: State) -> State:
    """Reflective y boundaries and periodic x boundaries."""
    s.rho = np.maximum(s.rho, 0.0)
    for arr in (s.rho, s.vx, s.vz, s.bz):
        arr[:, 0] = arr[:, 1]
        arr[:, -1] = arr[:, -2]
    s.vy[:, 0] = 0.0
    s.vy[:, -1] = 0.0

    # Extrapolate Az rather than forcing a constant wall value. This preserves
    # the physical wall slope Bx=dAz/dy and avoids a numerical current sheet.
    s.az[:, 0] = 2.0 * s.az[:, 1] - s.az[:, 2]
    s.az[:, -1] = 2.0 * s.az[:, -2] - s.az[:, -3]
    return s


def initial_state(x: np.ndarray, y: np.ndarray, p: Params) -> State:
    X, Y = np.meshgrid(x, y, indexing="ij")
    kx = np.pi / p.lx
    ky = np.pi / p.ly
    rho = 0.2 + 1.0 / np.cosh(Y / p.width) ** 2
    az = p.b0 * p.width * np.log(np.cosh(Y / p.width))
    az += p.psi0 * np.cos(ky * Y) * np.cos(2.0 * kx * X)
    az -= np.mean(az)
    zeros = np.zeros_like(rho)
    return apply_boundaries(
        State(
            rho=rho,
            vx=zeros.copy(),
            vy=zeros.copy(),
            vz=zeros.copy(),
            az=az,
            bz=zeros.copy(),
        )
    )


def magnetic_field(s: State, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bx = ddy(s.az, dy)
    by = -ddx(s.az, dx)
    return bx, by, s.bz


def current_density(
    bx: np.ndarray, by: np.ndarray, bz: np.ndarray, dx: float, dy: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    jx = ddy(bz, dy)
    jy = -ddx(bz, dx)
    jz = ddx(by, dx) - ddy(bx, dy)
    return jx, jy, jz


def rhs(s: State, p: Params, dx: float, dy: float) -> State:
    rho = np.maximum(s.rho, p.rho_floor)
    bx, by, bz = magnetic_field(s, dx, dy)
    jx, jy, jz = current_density(bx, by, bz, dx, dy)

    rhovx = rho * s.vx
    rhovy = rho * s.vy
    drho = -(ddx(rhovx, dx) + ddy(rhovy, dy)) + 0.25 * p.nu * lap(s.rho, dx, dy) - p.nu_h * hyp(s.rho, dx, dy)

    adv_vx = s.vx * ddx(s.vx, dx) + s.vy * ddy(s.vx, dy)
    adv_vy = s.vx * ddx(s.vy, dx) + s.vy * ddy(s.vy, dy)
    adv_vz = s.vx * ddx(s.vz, dx) + s.vy * ddy(s.vz, dy)

    lorentz_x = jy * bz - jz * by
    lorentz_y = jz * bx - jx * bz
    lorentz_z = jx * by - jy * bx

    dvx = -adv_vx - p.cs2 * ddx(rho, dx) / rho + lorentz_x / rho + p.nu * lap(s.vx, dx, dy) - p.nu_h * hyp(s.vx, dx, dy)
    dvy = -adv_vy - p.cs2 * ddy(rho, dy) / rho + lorentz_y / rho + p.nu * lap(s.vy, dx, dy) - p.nu_h * hyp(s.vy, dx, dy)
    dvz = -adv_vz + lorentz_z / rho + p.nu * lap(s.vz, dx, dy) - p.nu_h * hyp(s.vz, dx, dy)

    vx_bx = s.vy * bz - s.vz * by
    vx_by = s.vz * bx - s.vx * bz
    vx_bz = s.vx * by - s.vy * bx

    jxb_x = jy * bz - jz * by
    jxb_y = jz * bx - jx * bz
    jxb_z = jx * by - jy * bx

    fx = vx_bx - p.hall_coeff * jxb_x / rho
    fy = vx_by - p.hall_coeff * jxb_y / rho
    fz = vx_bz - p.hall_coeff * jxb_z / rho

    daz = fz + p.eta * lap(s.az, dx, dy) - p.eta_h * hyp(s.az, dx, dy)
    dbz = ddx(fy, dx) - ddy(fx, dy) + p.eta * lap(s.bz, dx, dy) - p.eta_h * hyp(s.bz, dx, dy)

    return State(drho, dvx, dvy, dvz, daz, dbz)


def add_scaled(s: State, k: State, scale: float, p: Params) -> State:
    out = State(
        rho=np.maximum(s.rho + scale * k.rho, p.rho_floor),
        vx=s.vx + scale * k.vx,
        vy=s.vy + scale * k.vy,
        vz=s.vz + scale * k.vz,
        az=s.az + scale * k.az,
        bz=s.bz + scale * k.bz,
    )
    return apply_boundaries(out)


def rk2_step(s: State, dt: float, p: Params, dx: float, dy: float) -> State:
    k1 = rhs(s, p, dx, dy)
    predictor = add_scaled(s, k1, dt, p)
    k2 = rhs(predictor, p, dx, dy)
    out = State(
        rho=np.maximum(0.5 * (s.rho + predictor.rho + dt * k2.rho), p.rho_floor),
        vx=0.5 * (s.vx + predictor.vx + dt * k2.vx),
        vy=0.5 * (s.vy + predictor.vy + dt * k2.vy),
        vz=0.5 * (s.vz + predictor.vz + dt * k2.vz),
        az=0.5 * (s.az + predictor.az + dt * k2.az),
        bz=0.5 * (s.bz + predictor.bz + dt * k2.bz),
    )
    return apply_boundaries(out)


def timestep(s: State, p: Params, dx: float, dy: float, remaining: float) -> float:
    bx, by, bz = magnetic_field(s, dx, dy)
    rho = np.maximum(s.rho, p.rho_floor)
    b2 = bx * bx + by * by + bz * bz
    vmag = np.sqrt(s.vx * s.vx + s.vy * s.vy + s.vz * s.vz)
    fast = np.sqrt(p.cs2 + b2 / rho) + vmag
    dt_mhd = p.cfl * min(dx, dy) / max(float(np.max(fast)), 1.0e-12)
    whistler = p.hall_coeff * np.sqrt(b2) / rho
    dt_hall = p.cfl_hall * min(dx, dy) ** 2 / max(float(np.max(whistler)), 1.0e-12)
    dt_diff = 0.2 * min(dx, dy) ** 2 / max(p.eta + p.nu, 1.0e-12)
    dt_hyper = 0.2 * min(dx, dy) ** 4 / max(p.eta_h + p.nu_h, 1.0e-12)
    return min(dt_mhd, dt_hall, dt_diff, dt_hyper, remaining)


def div_b(bx: np.ndarray, by: np.ndarray, dx: float, dy: float) -> np.ndarray:
    return ddx(bx, dx) + ddy(by, dy)


def force_balance_residual(s: State, p: Params, dx: float, dy: float) -> np.ndarray:
    bx, by, bz = magnetic_field(s, dx, dy)
    jx, _, jz = current_density(bx, by, bz, dx, dy)
    pressure_force_y = -p.cs2 * ddy(np.maximum(s.rho, p.rho_floor), dy)
    lorentz_force_y = jz * bx - jx * bz
    return pressure_force_y + lorentz_force_y


def reconnected_flux(by: np.ndarray, x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    j0 = int(np.argmin(np.abs(y)))
    xmax = 0.5 * (x.max() - x.min())
    mask = (x >= 0.0) & (x <= xmax)
    signed = float(_trapezoid(by[mask, j0], x[mask]))
    unsigned = float(_trapezoid(np.abs(by[mask, j0]), x[mask]))
    return signed, unsigned


def diagnostics(s: State, p: Params, x: np.ndarray, y: np.ndarray, dx: float, dy: float, t: float, step: int) -> dict[str, float]:
    bx, by, bz = magnetic_field(s, dx, dy)
    jx, jy, jz = current_density(bx, by, bz, dx, dy)
    signed, unsigned = reconnected_flux(by, x, y)
    bmag = np.sqrt(bx * bx + by * by + bz * bz)
    div = div_b(bx, by, dx, dy)
    force_res = force_balance_residual(s, p, dx, dy)
    rho = np.maximum(s.rho, 1.0e-12)
    hall_z = (jx * by - jy * bx) / rho
    return {
        "time": float(t),
        "step": int(step),
        "flux_signed": signed,
        "flux_unsigned": unsigned,
        "jz_abs_max": float(np.max(np.abs(jz))),
        "jz_l2": float(np.sqrt(np.mean(jz * jz))),
        "divB_l2": float(np.sqrt(np.mean(div * div))),
        "force_balance_l2": float(np.sqrt(np.mean(force_res * force_res))),
        "force_balance_abs_max": float(np.max(np.abs(force_res))),
        "rho_min": float(np.min(s.rho)),
        "rho_max": float(np.max(s.rho)),
        "B_abs_max": float(np.max(bmag)),
        "Bz_abs_max": float(np.max(np.abs(bz))),
        "vz_abs_max": float(np.max(np.abs(s.vz))),
        "hall_z_abs_max": float(np.max(np.abs(hall_z))),
        "kinetic_energy": float(0.5 * np.mean(rho * (s.vx * s.vx + s.vy * s.vy + s.vz * s.vz))),
        "magnetic_energy": float(0.5 * np.mean(bmag * bmag)),
    }


def save_snapshot(path: Path, s: State, x: np.ndarray, y: np.ndarray, dx: float, dy: float, t: float, step: int) -> None:
    bx, by, bz = magnetic_field(s, dx, dy)
    np.savez_compressed(
        path,
        time=t,
        step=step,
        x=x,
        y=y,
        rho=s.rho,
        vx=s.vx,
        vy=s.vy,
        vz=s.vz,
        Bx=bx,
        By=by,
        Bz=bz,
        Az=s.az,
    )


def plot_state(path: Path, s: State, x: np.ndarray, y: np.ndarray, dx: float, dy: float, t: float) -> None:
    bx, by, bz = magnetic_field(s, dx, dy)
    _, _, jz = current_density(bx, by, bz, dx, dy)
    div = div_b(bx, by, dx, dy)
    fields = [
        (s.rho, r"$\rho$", "viridis", False),
        (s.vx, r"$v_x$", "RdBu_r", True),
        (s.vy, r"$v_y$", "RdBu_r", True),
        (bx, r"$B_x$", "RdBu_r", True),
        (by, r"$B_y$", "RdBu_r", True),
        (bz, r"$B_z$", "RdBu_r", True),
        (jz, r"$J_z$", "RdBu_r", True),
        (np.sqrt(bx * bx + by * by + bz * bz), r"$|B|$", "magma", False),
        (div, r"$\nabla\cdot B$", "RdBu_r", True),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(18, 11))
    fig.suptitle(f"Python Hall-MHD Harris current sheet (t={t:.3f})")
    for ax, (data, label, cmap, symmetric) in zip(axes.flat, fields):
        norm = None
        vmin = vmax = None
        if label in CLIMS:
            vmin, vmax = CLIMS[label]
        elif symmetric:
            vmax = max(float(np.percentile(np.abs(data), 99.5)), 1.0e-12)
            vmin = -vmax
            if vmax / max(float(np.max(np.abs(data))), 1.0e-12) < 0.9:
                norm = SymLogNorm(linthresh=vmax / 100.0, vmin=vmin, vmax=vmax)
                vmin = vmax = None
        mesh = ax.pcolormesh(x, y, data.T, shading="auto", cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
        fig.colorbar(mesh, ax=ax, shrink=0.82)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(label)
        ax.set_aspect("equal")
    for ax in axes.flat[len(fields):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_timeseries(path: Path, rows: list[dict[str, float]]) -> None:
    t = np.array([r["time"] for r in rows])
    flux = np.array([r["flux_unsigned"] for r in rows])
    jmax = np.array([r["jz_abs_max"] for r in rows])
    div = np.array([r["divB_l2"] for r in rows])
    bz = np.array([r["Bz_abs_max"] for r in rows])
    force = np.array([r["force_balance_l2"] for r in rows])
    rho_min = np.array([r["rho_min"] for r in rows])

    fig, axes = plt.subplots(3, 2, figsize=(12, 11))
    series = [
        (flux, r"$\int_0^{L_x/2}|B_y(x,0)|dx$", "C0"),
        (jmax, r"$\max |J_z|$", "C3"),
        (div, r"$||\nabla\cdot B||_2$", "C4"),
        (bz, r"$\max |B_z|$", "C2"),
        (force, r"$||-\nabla p + J\times B||_2$", "C1"),
        (rho_min, r"$\min \rho$", "C5"),
    ]
    for idx, (ax, (values, label, color)) in enumerate(zip(axes.flat, series)):
        if idx == 2 and np.all(values > 0.0):
            ax.semilogy(t, values, "o-", color=color)
        else:
            ax.plot(t, values, "o-", color=color)
        ax.set_xlabel("t")
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def make_gif(output_dir: Path) -> Path:
    """Create an animated GIF from saved .npz snapshots."""
    from PIL import Image

    snapshots = sorted(output_dir.glob("python_hall_mhd_*.npz"))
    if not snapshots:
        raise FileNotFoundError("No .npz snapshots found")

    frames: list[Image.Image] = []
    for snp_path in snapshots:
        data = np.load(snp_path)
        bx, by = data["Bx"], data["By"]
        jz = np.gradient(by, data["x"], axis=0, edge_order=2) - np.gradient(bx, data["y"], axis=1, edge_order=2)
        t = float(data["time"])

        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        for ax, (arr, label, cmap, sym) in zip(
            axes,
            [
                (data["rho"], r"$\rho$", "viridis", False),
                (jz, r"$J_z$", "RdBu_r", True),
                (np.sqrt(bx * bx + by * by + data["Bz"] ** 2), r"$|\mathbf{B}|$", "magma", False),
            ],
        ):
            if sym:
                vmax = max(float(np.percentile(np.abs(arr), 99.5)), 1e-12)
                mesh = ax.pcolormesh(data["x"], data["y"], arr.T, shading="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
            else:
                mesh = ax.pcolormesh(data["x"], data["y"], arr.T, shading="auto", cmap=cmap)
            fig.colorbar(mesh, ax=ax, shrink=0.85)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_title(f"{label}  t={t:.1f}")
            ax.set_aspect("equal")
        fig.suptitle("Python Hall-MHD Harris sheet", fontsize=13)
        fig.tight_layout()

        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())
        frames.append(Image.fromarray(buf))
        plt.close(fig)

    gif_path = output_dir / "evolution.gif"
    frames[0].save(gif_path, save_all=True, append_images=frames[1:], loop=0, duration=600)
    return gif_path


def write_diagnostics(path: Path, rows: list[dict[str, float]]) -> None:
    columns = [
        "time",
        "step",
        "flux_signed",
        "flux_unsigned",
        "jz_abs_max",
        "jz_l2",
        "divB_l2",
        "force_balance_l2",
        "force_balance_abs_max",
        "rho_min",
        "rho_max",
        "B_abs_max",
        "Bz_abs_max",
        "vz_abs_max",
        "hall_z_abs_max",
        "kinetic_energy",
        "magnetic_energy",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_metadata(path: Path, p: Params, dx: float, dy: float, rows: list[dict[str, float]]) -> None:
    params = asdict(p)
    params["output_dir"] = str(p.output_dir)
    final = rows[-1] if rows else {}
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv,
        "parameters": params,
        "grid": {
            "dx": dx,
            "dy": dy,
            "cells_per_sheet_width": p.width / dy,
            "domain": {"lx": p.lx, "ly": p.ly},
        },
        "outputs": {
            "snapshot_count": len(rows),
            "diagnostics_csv": "python_hall_mhd_diagnostics.csv",
            "timeseries_png": "python_hall_mhd_timeseries.png",
            "gif": "evolution.gif",
        },
        "final_diagnostics": final,
        "notes": [
            "Python finite-difference Hall-MHD run; PLUTO remains the quantitative reference.",
            "The explicit Hall timestep scales approximately as dx^2.",
        ],
    }
    path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def run_simulation(p: Params) -> list[dict[str, float]]:
    out = p.output_dir
    out.mkdir(parents=True, exist_ok=True)
    for child in out.glob("python_hall_mhd_*"):
        if child.is_file():
            child.unlink()

    x, y, dx, dy = make_grid(p)
    state = initial_state(x, y, p)
    initial_force = force_balance_residual(state, p, dx, dy)
    print(
        "Initial force balance residual: "
        f"L2={np.sqrt(np.mean(initial_force * initial_force)):.6e}, "
        f"max={np.max(np.abs(initial_force)):.6e}"
    )
    rows: list[dict[str, float]] = []
    t = 0.0
    step = 0
    next_output = p.output_dt
    snap_id = 0

    def emit() -> None:
        nonlocal snap_id
        rows.append(diagnostics(state, p, x, y, dx, dy, t, step))
        save_snapshot(out / f"python_hall_mhd_{snap_id:04d}.npz", state, x, y, dx, dy, t, step)
        plot_state(out / f"python_hall_mhd_{snap_id:04d}_t{t:.2f}.png", state, x, y, dx, dy, t)
        print(f"Snapshot {snap_id:04d}: t={t:.6f}, step={step}")
        snap_id += 1

    emit()
    while t < p.tstop - 1.0e-12:
        if step >= p.max_steps:
            raise RuntimeError(f"Reached max_steps={p.max_steps} before tstop={p.tstop}")
        target = min(p.tstop, next_output if next_output > t + 1.0e-12 else p.tstop)
        dt = timestep(state, p, dx, dy, target - t)
        state = rk2_step(state, dt, p, dx, dy)
        if not all(np.all(np.isfinite(arr)) for arr in (state.rho, state.vx, state.vy, state.vz, state.az, state.bz)):
            raise FloatingPointError(f"Non-finite values at t={t:.6f}, step={step}")
        t += dt
        step += 1
        if step % 5000 == 0:
            print(f"Progress: t={t:.6f}/{p.tstop:.6f}, step={step}, dt={dt:.3e}")
        if t >= next_output - 1.0e-10:
            emit()
            next_output += p.output_dt

    if rows[-1]["time"] < p.tstop - 1.0e-10:
        emit()
    write_diagnostics(out / "python_hall_mhd_diagnostics.csv", rows)
    plot_timeseries(out / "python_hall_mhd_timeseries.png", rows)
    write_metadata(out / "run_metadata.json", p, dx, dy, rows)
    try:
        gif = make_gif(out)
        print(f"Animation saved: {gif}")
    except (ImportError, FileNotFoundError) as e:
        print(f"GIF not generated ({e})")
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=64, help="Number of x cells.")
    parser.add_argument("--ny", type=int, default=32, help="Number of y cells.")
    parser.add_argument("--tstop", type=float, default=5.0, help="Final time.")
    parser.add_argument("--output-dt", type=float, default=1.0, help="Time between saved snapshots.")
    parser.add_argument("--hall-coeff", type=float, default=1.0, help="Hall coefficient d_i in normalized units.")
    parser.add_argument("--cs2", type=float, default=0.5, help="Isothermal sound speed squared.")
    parser.add_argument("--b0", type=float, default=1.0, help="Asymptotic Harris-sheet magnetic field.")
    parser.add_argument("--width", type=float, default=0.5, help="Harris current-sheet half-width.")
    parser.add_argument("--psi0", type=float, default=0.02, help="Initial magnetic perturbation amplitude.")
    parser.add_argument("--eta", type=float, default=5.0e-3, help="Magnetic diffusivity used for numerical regularization.")
    parser.add_argument("--nu", type=float, default=1.0e-3, help="Kinematic viscosity used for numerical regularization.")
    parser.add_argument("--eta-h", type=float, default=1.0e-4, help="Hyper-resistivity (4th order) for grid-scale damping.")
    parser.add_argument("--nu-h", type=float, default=5.0e-5, help="Hyper-viscosity (4th order) for grid-scale damping.")
    parser.add_argument("--cfl", type=float, default=0.22, help="MHD CFL factor.")
    parser.add_argument("--cfl-hall", type=float, default=0.18, help="Hall/whistler CFL factor.")
    parser.add_argument("--max-steps", type=int, default=1_000_000, help="Safety limit for explicit time steps.")
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR, help="Directory for figures, snapshots and diagnostics.")
    parser.add_argument(
        "--pluto-grid",
        action="store_true",
        help="Use the PLUTO grid and end time: nx=256, ny=128, tstop=60, output_dt=5.",
    )
    parser.add_argument(
        "--long-run",
        action="store_true",
        help="Use conservative damping for longer exploratory runs: tstop=15, output_dt=3.",
    )
    return parser.parse_args()


def params_from_args(args: argparse.Namespace) -> Params:
    p = Params(
        nx=args.nx,
        ny=args.ny,
        tstop=args.tstop,
        output_dt=args.output_dt,
        cs2=args.cs2,
        b0=args.b0,
        width=args.width,
        hall_coeff=args.hall_coeff,
        psi0=args.psi0,
        eta=args.eta,
        nu=args.nu,
        eta_h=args.eta_h,
        nu_h=args.nu_h,
        cfl=args.cfl,
        cfl_hall=args.cfl_hall,
        max_steps=args.max_steps,
        output_dir=args.output_dir,
    )
    if args.pluto_grid:
        p = replace(p, nx=256, ny=128, tstop=60.0, output_dt=5.0)
    if args.long_run:
        p = replace(p, tstop=15.0, output_dt=3.0, eta=2.0e-2, nu=1.0e-2, eta_h=1.0e-4, nu_h=5.0e-5)
    return p


def main() -> None:
    args = parse_args()
    p = params_from_args(args)
    print(
        "Running Python Hall-MHD Harris sheet: "
        f"{p.nx}x{p.ny}, tstop={p.tstop}, output_dt={p.output_dt}, hall={p.hall_coeff},"
        f" psi0={p.psi0}, eta={p.eta}, nu={p.nu}, eta_h={p.eta_h}, nu_h={p.nu_h}"
    )
    rows = run_simulation(p)
    final = rows[-1]
    print(f"Wrote Python simulation products to {p.output_dir}")
    print(f"Final time: {final['time']:.6f} at step {final['step']}")
    print(f"Final unsigned reconnected flux: {final['flux_unsigned']:.6f}")
    print(f"Final max |Jz|: {final['jz_abs_max']:.6f}")
    print(f"Final max |Bz|: {final['Bz_abs_max']:.6f}")
    print(f"Final ||div B||_2: {final['divB_l2']:.6e}")


if __name__ == "__main__":
    main()
