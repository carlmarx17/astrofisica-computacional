"""
Post-processing for the Hall-MHD Harris current sheet project.

This script loads the PLUTO VTK snapshots, rebuilds the same Harris-sheet
initial condition in Python, computes quantitative diagnostics, and writes the
figures used in the report.
"""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / 'figs'
PLUTO_RUN = ROOT / 'pluto_sim'
OUT.mkdir(parents=True, exist_ok=True)

pluto_py_candidates = []
if os.environ.get('PLUTO_DIR'):
    pluto_py_candidates.append(Path(os.environ['PLUTO_DIR']) / 'Tools' / 'pyPLUTO')
pluto_py_candidates.append(ROOT.parent / 'PLUTO' / 'Tools' / 'pyPLUTO')
for candidate in pluto_py_candidates:
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))
        break
import pyPLUTO.pload as pp

CS2 = 0.5
B0 = 1.0
WIDTH = 0.5
PSI0 = 0.02


def read_snapshots(run_dir):
    snapshots = []
    with open(run_dir / 'vtk.out') as fh:
        for line in fh:
            parts = line.split()
            snapshots.append({'n': int(parts[0]), 'time': float(parts[1]), 'step': int(parts[3])})
    return snapshots


def load_pluto_snapshot(run_dir, snapshot):
    d = pp.pload(snapshot['n'], w_dir=str(run_dir) + '/', datatype='vtk')
    return {
        'n': snapshot['n'],
        'time': snapshot['time'],
        'step': snapshot['step'],
        'x': np.asarray(d.x1),
        'y': np.asarray(d.x2),
        'rho': np.asarray(d.rho),
        'vx': np.asarray(d.vx1),
        'vy': np.asarray(d.vx2),
        'Bx': np.asarray(d.Bx1),
        'By': np.asarray(d.Bx2),
    }


def harris_initial(x, y):
    X, Y = np.meshgrid(x, y, indexing='ij')
    lx = x.max() - x.min() + (x[1] - x[0])
    ly = y.max() - y.min() + (y[1] - y[0])
    kx = np.pi / lx
    ky = np.pi / ly

    rho = 0.2 + 1.0 / np.cosh(Y / WIDTH) ** 2
    Bx = B0 * np.tanh(Y / WIDTH)
    By = np.zeros_like(Bx)
    Bx += -PSI0 * ky * np.sin(ky * Y) * np.cos(2.0 * kx * X)
    By += PSI0 * 2.0 * kx * np.sin(2.0 * kx * X) * np.cos(ky * Y)
    return {
        'rho': rho,
        'vx': np.zeros_like(rho),
        'vy': np.zeros_like(rho),
        'Bx': Bx,
        'By': By,
    }


def grid_spacing(values):
    return float(np.mean(np.diff(values)))


def current_jz(Bx, By, x, y):
    dby_dx = np.gradient(By, x, axis=0, edge_order=2)
    dbx_dy = np.gradient(Bx, y, axis=1, edge_order=2)
    return dby_dx - dbx_dy


def div_b(Bx, By, x, y):
    dbx_dx = np.gradient(Bx, x, axis=0, edge_order=2)
    dby_dy = np.gradient(By, y, axis=1, edge_order=2)
    return dbx_dx + dby_dy


def reconnected_flux_midplane(By, x, y):
    j0 = int(np.argmin(np.abs(y)))
    xmax = 0.5 * (x.max() - x.min())
    mask = (x >= 0.0) & (x <= xmax)
    signed = np.trapezoid(By[mask, j0], x[mask])
    unsigned = np.trapezoid(np.abs(By[mask, j0]), x[mask])
    return signed, unsigned


def rel_l2(a, b):
    denom = np.sqrt(np.mean(np.asarray(b) ** 2))
    if denom < 1e-14:
        return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)) / denom)


def write_csv(path, rows, columns):
    with open(path, 'w') as fh:
        fh.write(','.join(columns) + '\n')
        for row in rows:
            fh.write(','.join(f'{row[col]:.10e}' if isinstance(row[col], float) else str(row[col]) for col in columns) + '\n')


def plot_timeseries(rows):
    t = np.array([r['time'] for r in rows])
    flux = np.array([r['flux_unsigned'] for r in rows])
    jmax = np.array([r['jz_abs_max'] for r in rows])
    div_l2 = np.array([r['divB_l2'] for r in rows])

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    axes[0].plot(t, flux, 'o-', color='C0')
    axes[0].set_xlabel('t')
    axes[0].set_ylabel(r'$\int_0^{L_x/2} |B_y(x,0)| dx$')
    axes[0].set_title('Reconnected flux')

    axes[1].plot(t, jmax, 's-', color='C3')
    axes[1].set_xlabel('t')
    axes[1].set_ylabel(r'$\max |J_z|$')
    axes[1].set_title('Current sheet intensity')

    axes[2].semilogy(t, div_l2, '^-', color='C4')
    axes[2].set_xlabel('t')
    axes[2].set_ylabel(r'$||\nabla\cdot B||_2$')
    axes[2].set_title('Divergence error')

    for ax in axes:
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / 'diagnostics_timeseries.png', dpi=140)
    plt.close(fig)


def plot_final_state(d):
    fields = [
        ('rho', d['rho'], r'$\rho$', 'plasma', False),
        ('pressure', CS2 * d['rho'], r'$P=c_s^2\rho$', 'inferno', False),
        ('vx', d['vx'], r'$v_x$', 'RdBu_r', True),
        ('vy', d['vy'], r'$v_y$', 'RdBu_r', True),
        ('Bx', d['Bx'], r'$B_x$', 'RdBu_r', True),
        ('By', d['By'], r'$B_y$', 'RdBu_r', True),
        ('Bmag', np.sqrt(d['Bx'] ** 2 + d['By'] ** 2), r'$|B|$', 'magma', False),
        ('Jz', current_jz(d['Bx'], d['By'], d['x'], d['y']), r'$J_z$', 'RdBu_r', True),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    fig.suptitle(f'PLUTO Hall MHD Harris Current Sheet (t={d["time"]:.1f})', fontsize=16)
    for ax, (_, data, label, cmap, symmetric) in zip(axes.flat, fields):
        norm = None
        if symmetric and np.max(np.abs(data)) > 1e-14:
            vmax = float(np.max(np.abs(data)))
            norm = SymLogNorm(linthresh=vmax / 100.0, vmin=-vmax, vmax=vmax)
        mesh = ax.pcolormesh(d['x'], d['y'], data.T, shading='auto', cmap=cmap, norm=norm)
        fig.colorbar(mesh, ax=ax, shrink=0.8)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(label)
        ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(OUT / 'pluto_final_all_variables.png', dpi=140)
    plt.close(fig)


def plot_jz_fieldlines(d):
    jz = current_jz(d['Bx'], d['By'], d['x'], d['y'])
    vmax = float(np.percentile(np.abs(jz), 99.5))
    fig, ax = plt.subplots(figsize=(10, 5))
    mesh = ax.pcolormesh(d['x'], d['y'], jz.T, shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    fig.colorbar(mesh, ax=ax, label=r'$J_z = \partial_x B_y - \partial_y B_x$')
    ax.streamplot(d['x'], d['y'], d['Bx'].T, d['By'].T, color='k', density=1.4, linewidth=0.65, arrowsize=0.8)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Current density and magnetic field lines (nearest to t=57: t={d["time"]:.1f})')
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(OUT / f'jz_fieldlines_t{d["time"]:.0f}.png', dpi=160)
    plt.close(fig)


def plot_python_pluto_initial(d0, analytic):
    variables = [
        ('rho', r'$\rho$'),
        ('Bx', r'$B_x$'),
        ('By', r'$B_y$'),
    ]
    fig, axes = plt.subplots(len(variables), 3, figsize=(15, 11))
    for row, (name, label) in enumerate(variables):
        data = [analytic[name], d0[name], d0[name] - analytic[name]]
        titles = [f'Python initial {label}', f'PLUTO t=0 {label}', 'PLUTO - Python']
        for col, (arr, title) in enumerate(zip(data, titles)):
            ax = axes[row, col]
            cmap = 'RdBu_r' if col == 2 or np.nanmin(arr) < 0 else 'viridis'
            vmax = None
            vmin = None
            if col == 2:
                vmax = max(float(np.max(np.abs(arr))), 1e-14)
                vmin = -vmax
            mesh = ax.pcolormesh(d0['x'], d0['y'], arr.T, shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
            fig.colorbar(mesh, ax=ax, shrink=0.8)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_title(title)
            ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(OUT / 'python_pluto_initial_comparison.png', dpi=140)
    plt.close(fig)


def main():
    snapshots = read_snapshots(PLUTO_RUN)
    data = [load_pluto_snapshot(PLUTO_RUN, s) for s in snapshots]
    analytic0 = harris_initial(data[0]['x'], data[0]['y'])

    rows = []
    for d in data:
        jz = current_jz(d['Bx'], d['By'], d['x'], d['y'])
        divergence = div_b(d['Bx'], d['By'], d['x'], d['y'])
        signed, unsigned = reconnected_flux_midplane(d['By'], d['x'], d['y'])
        bmag = np.sqrt(d['Bx'] ** 2 + d['By'] ** 2)
        rows.append({
            'snapshot': d['n'],
            'time': float(d['time']),
            'step': d['step'],
            'flux_signed': float(signed),
            'flux_unsigned': float(unsigned),
            'jz_abs_max': float(np.max(np.abs(jz))),
            'jz_l2': float(np.sqrt(np.mean(jz ** 2))),
            'divB_l2': float(np.sqrt(np.mean(divergence ** 2))),
            'rho_min': float(np.min(d['rho'])),
            'rho_max': float(np.max(d['rho'])),
            'B_abs_max': float(np.max(bmag)),
        })

    initial_error_rows = []
    for name in ['rho', 'vx', 'vy', 'Bx', 'By']:
        initial_error_rows.append({
            'variable': name,
            'relative_l2_error': rel_l2(data[0][name], analytic0[name]),
            'max_abs_error': float(np.max(np.abs(data[0][name] - analytic0[name]))),
        })

    write_csv(
        OUT / 'diagnostics.csv',
        rows,
        ['snapshot', 'time', 'step', 'flux_signed', 'flux_unsigned', 'jz_abs_max', 'jz_l2', 'divB_l2', 'rho_min', 'rho_max', 'B_abs_max'],
    )
    write_csv(
        OUT / 'initial_condition_errors.csv',
        initial_error_rows,
        ['variable', 'relative_l2_error', 'max_abs_error'],
    )

    nearest_57 = min(data, key=lambda d: abs(d['time'] - 57.0))
    plot_timeseries(rows)
    plot_final_state(data[-1])
    plot_jz_fieldlines(nearest_57)
    plot_python_pluto_initial(data[0], analytic0)

    print('Wrote analysis products to', OUT)
    print('Final unsigned reconnected flux:', f'{rows[-1]["flux_unsigned"]:.6f}')
    print('Nearest snapshot to t=57:', nearest_57['n'], f't={nearest_57["time"]:.3f}')


if __name__ == '__main__':
    main()
