"""Generate per-snapshot PLUTO panels for the Harris current sheet."""
from pathlib import Path
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / 'pluto_sim'
PLOT_DIR = ROOT / 'plots'
PLOT_DIR.mkdir(parents=True, exist_ok=True)

candidates = []
if os.environ.get('PLUTO_DIR'):
    candidates.append(Path(os.environ['PLUTO_DIR']) / 'Tools' / 'pyPLUTO')
candidates.append(ROOT.parent / 'PLUTO' / 'Tools' / 'pyPLUTO')
for candidate in candidates:
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))
        break
import pyPLUTO.pload as pp

CS2 = 0.5


def current_jz(Bx, By, x, y):
    return np.gradient(By, x, axis=0, edge_order=2) - np.gradient(Bx, y, axis=1, edge_order=2)


def symmetric_norm(data):
    if np.nanmin(data) < 0 and np.nanmax(data) > 0 and np.nanmax(np.abs(data)) > 1e-14:
        vmax = float(np.nanmax(np.abs(data)))
        return SymLogNorm(linthresh=vmax / 100.0, vmin=-vmax, vmax=vmax)
    return None


def read_snapshots():
    snapshots = []
    with open(RUN_DIR / 'vtk.out') as fh:
        for line in fh:
            parts = line.split()
            snapshots.append({'n': int(parts[0]), 'time': float(parts[1])})
    return snapshots


def main():
    snapshots = read_snapshots()
    print(f'Found {len(snapshots)} VTK snapshots')

    for snapshot in snapshots:
        n = snapshot['n']
        t = snapshot['time']
        print(f'  Plotting snapshot {n}: t = {t:.2f}')
        d = pp.pload(n, w_dir=str(RUN_DIR) + '/', datatype='vtk')

        Bmag = np.sqrt(d.Bx1 ** 2 + d.Bx2 ** 2)
        Jz = current_jz(d.Bx1, d.Bx2, d.x1, d.x2)
        fields = [
            (d.rho, r'$\rho$', 'plasma', None),
            (CS2 * d.rho, r'$P=c_s^2\rho$', 'inferno', None),
            (d.vx1, r'$v_x$', 'RdBu_r', symmetric_norm(d.vx1)),
            (d.vx2, r'$v_y$', 'RdBu_r', symmetric_norm(d.vx2)),
            (d.Bx1, r'$B_x$', 'RdBu_r', symmetric_norm(d.Bx1)),
            (d.Bx2, r'$B_y$', 'RdBu_r', symmetric_norm(d.Bx2)),
            (Bmag, r'$|B|$', 'magma', None),
            (Jz, r'$J_z$', 'RdBu_r', symmetric_norm(Jz)),
        ]

        fig, axes = plt.subplots(2, 4, figsize=(20, 9))
        fig.suptitle(f'Hall MHD - Harris Current Sheet (t = {t:.1f})', fontsize=16)
        for ax, (data, label, cmap, norm) in zip(axes.flat, fields):
            im = ax.pcolormesh(d.x1, d.x2, data.T, shading='auto', norm=norm, cmap=cmap)
            fig.colorbar(im, ax=ax, shrink=0.8)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_title(label)
            ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f'hall_cs_{n:04d}_t{t:.1f}.png', dpi=120)
        plt.close(fig)

    print('All plots saved in', PLOT_DIR)


if __name__ == '__main__':
    main()
