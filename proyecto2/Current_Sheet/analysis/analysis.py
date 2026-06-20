"""
Post-processing for the Hall-MHD Harris current sheet project.

This script loads the PLUTO VTK snapshots, rebuilds the same Harris-sheet
initial condition in Python, computes quantitative diagnostics, and writes the
figures used in the report.
"""
from pathlib import Path
import os
import sys

MPLCONFIG = Path('/tmp') / 'matplotlib-cache'
MPLCONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(MPLCONFIG))

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
        'vz': np.asarray(d.vx3),
        'Bx': np.asarray(d.Bx1),
        'By': np.asarray(d.Bx2),
        'Bz': np.asarray(d.Bx3),
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
        'vz': np.zeros_like(rho),
        'Bx': Bx,
        'By': By,
        'Bz': np.zeros_like(rho),
    }


def grid_spacing(values):
    return float(np.mean(np.diff(values)))


def current_jz(Bx, By, x, y):
    dby_dx = np.gradient(By, x, axis=0, edge_order=2)
    dbx_dy = np.gradient(Bx, y, axis=1, edge_order=2)
    return dby_dx - dbx_dy


def current_components(Bx, By, Bz, x, y):
    jx = np.gradient(Bz, y, axis=1, edge_order=2)
    jy = -np.gradient(Bz, x, axis=0, edge_order=2)
    jz = current_jz(Bx, By, x, y)
    return jx, jy, jz


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


def vector_potential_az(Bx, By, x, y):
    """Reconstruct Az from Bx=dAz/dy and By=-dAz/dx using a fixed path."""
    az = np.zeros_like(Bx)
    for j in range(1, len(y)):
        dy = y[j] - y[j - 1]
        az[0, j] = az[0, j - 1] + 0.5 * (Bx[0, j] + Bx[0, j - 1]) * dy
    for j in range(len(y)):
        for i in range(1, len(x)):
            dx = x[i] - x[i - 1]
            az[i, j] = az[i - 1, j] - 0.5 * (By[i, j] + By[i - 1, j]) * dx
    return az


def ox_points_from_az(Az, x, y):
    j0 = int(np.argmin(np.abs(y)))
    i_o = int(np.argmin(np.abs(x)))
    mid = Az[:, j0]
    left = np.where((x > -0.45 * (x.max() - x.min())) & (x < -1.0))[0]
    right = np.where((x < 0.45 * (x.max() - x.min())) & (x > 1.0))[0]
    i_x_left = int(left[np.argmax(mid[left])]) if len(left) else i_o
    i_x_right = int(right[np.argmax(mid[right])]) if len(right) else i_o
    i_x = i_x_left if abs(Az[i_x_left, j0] - Az[i_o, j0]) > abs(Az[i_x_right, j0] - Az[i_o, j0]) else i_x_right
    return {
        'O': (i_o, j0, float(x[i_o]), float(y[j0]), float(Az[i_o, j0])),
        'X_left': (i_x_left, j0, float(x[i_x_left]), float(y[j0]), float(Az[i_x_left, j0])),
        'X_right': (i_x_right, j0, float(x[i_x_right]), float(y[j0]), float(Az[i_x_right, j0])),
        'X_used': (i_x, j0, float(x[i_x]), float(y[j0]), float(Az[i_x, j0])),
        'psi_az': float(abs(Az[i_o, j0] - Az[i_x, j0])),
    }


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
    az = vector_potential_az(d['Bx'], d['By'], d['x'], d['y'])
    points = ox_points_from_az(az, d['x'], d['y'])
    vmax = float(np.percentile(np.abs(jz), 99.5))
    fig, ax = plt.subplots(figsize=(10, 5))
    mesh = ax.pcolormesh(d['x'], d['y'], jz.T, shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    fig.colorbar(mesh, ax=ax, label=r'$J_z = \partial_x B_y - \partial_y B_x$')
    ax.streamplot(d['x'], d['y'], d['Bx'].T, d['By'].T, color='k', density=1.4, linewidth=0.65, arrowsize=0.8)
    for label, marker, color in [('O', 'o', 'gold'), ('X_left', 'x', 'lime'), ('X_right', 'x', 'lime')]:
        _, _, px, py, _ = points[label]
        if marker == 'x':
            ax.scatter(px, py, marker=marker, s=95, c=color, linewidths=1.4, zorder=5)
        else:
            ax.scatter(px, py, marker=marker, s=95, c=color, edgecolors='k', linewidths=0.8, zorder=5)
    ax.annotate('O-point', xy=(points['O'][2], points['O'][3]), xytext=(points['O'][2] + 0.8, points['O'][3] + 0.8),
                arrowprops={'arrowstyle': '->', 'color': 'k', 'lw': 1.0}, fontsize=9)
    ax.annotate('X-points', xy=(points['X_right'][2], points['X_right'][3]), xytext=(points['X_right'][2] - 0.6, points['X_right'][3] + 1.0),
                arrowprops={'arrowstyle': '->', 'color': 'k', 'lw': 1.0}, fontsize=9, ha='right')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Current density and field lines: O/X structure (t={d["time"]:.1f})')
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(OUT / f'jz_fieldlines_t{d["time"]:.0f}.png', dpi=160)
    plt.close(fig)


def plot_hall_signature(d):
    jx, jy, jz = current_components(d['Bx'], d['By'], d['Bz'], d['x'], d['y'])
    rho = np.maximum(d['rho'], 1e-12)
    hall_z = (jx * d['By'] - jy * d['Bx']) / rho
    electron_drift = np.sqrt((jx / rho) ** 2 + (jy / rho) ** 2 + (jz / rho) ** 2)

    fields = [
        (d['Bz'], r'$B_z$', 'RdBu_r'),
        (hall_z, r'$(J\times B)_z/\rho$', 'RdBu_r'),
        (electron_drift, r'$|\mathbf{v}_e-\mathbf{v}|=|\mathbf{J}|/\rho$', 'magma'),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, (data, label, cmap) in zip(axes, fields):
        kwargs = {}
        if np.nanmin(data) < 0 and np.nanmax(data) > 0:
            vmax = float(np.percentile(np.abs(data), 99.5))
            kwargs.update({'vmin': -vmax, 'vmax': vmax})
        mesh = ax.pcolormesh(d['x'], d['y'], data.T, shading='auto', cmap=cmap, **kwargs)
        fig.colorbar(mesh, ax=ax, shrink=0.82)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(label)
        ax.set_aspect('equal')
    fig.suptitle(f'Hall-specific diagnostics (t={d["time"]:.1f})')
    fig.tight_layout()
    fig.savefig(OUT / 'hall_signature_t60.png', dpi=150)
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
        jx, jy, jz = current_components(d['Bx'], d['By'], d['Bz'], d['x'], d['y'])
        divergence = div_b(d['Bx'], d['By'], d['x'], d['y'])
        signed, unsigned = reconnected_flux_midplane(d['By'], d['x'], d['y'])
        bmag = np.sqrt(d['Bx'] ** 2 + d['By'] ** 2 + d['Bz'] ** 2)
        rho = np.maximum(d['rho'], 1e-12)
        hall_z = (jx * d['By'] - jy * d['Bx']) / rho
        az = vector_potential_az(d['Bx'], d['By'], d['x'], d['y'])
        ox = ox_points_from_az(az, d['x'], d['y'])
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
            'Bz_abs_max': float(np.max(np.abs(d['Bz']))),
            'vz_abs_max': float(np.max(np.abs(d['vz']))),
            'hall_z_abs_max': float(np.max(np.abs(hall_z))),
            'psi_az': ox['psi_az'],
        })

    initial_error_rows = []
    for name in ['rho', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz']:
        initial_error_rows.append({
            'variable': name,
            'relative_l2_error': rel_l2(data[0][name], analytic0[name]),
            'max_abs_error': float(np.max(np.abs(data[0][name] - analytic0[name]))),
        })

    write_csv(
        OUT / 'diagnostics.csv',
        rows,
        ['snapshot', 'time', 'step', 'flux_signed', 'flux_unsigned', 'jz_abs_max', 'jz_l2', 'divB_l2', 'rho_min', 'rho_max', 'B_abs_max', 'Bz_abs_max', 'vz_abs_max', 'hall_z_abs_max', 'psi_az'],
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
    plot_hall_signature(data[-1])
    plot_python_pluto_initial(data[0], analytic0)

    az_final = vector_potential_az(data[-1]['Bx'], data[-1]['By'], data[-1]['x'], data[-1]['y'])
    points = ox_points_from_az(az_final, data[-1]['x'], data[-1]['y'])

    print('Wrote analysis products to', OUT)
    print('Final unsigned reconnected flux:', f'{rows[-1]["flux_unsigned"]:.6f}')
    print('Az O-X flux estimate:', f'{points["psi_az"]:.6f}')
    print('Nearest snapshot to t=57:', nearest_57['n'], f't={nearest_57["time"]:.3f}')


if __name__ == '__main__':
    main()
