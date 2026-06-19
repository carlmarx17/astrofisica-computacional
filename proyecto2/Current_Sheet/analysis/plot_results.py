import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

pluto_py = os.path.expandvars('$PLUTO_DIR/Tools/pyPLUTO')
sys.path.insert(0, pluto_py)
import pyPLUTO as pypl
import pyPLUTO.pload as pp

wdir = os.path.dirname(os.path.abspath(__file__)) + '/'

with open(wdir + 'vtk.out') as f:
    lines = f.readlines()

snapshots = []
for line in lines:
    parts = line.split()
    snapshots.append({'n': int(parts[0]), 'time': float(parts[1])})

print(f"Found {len(snapshots)} VTK snapshots")

for s in snapshots:
    n = s['n']
    t = s['time']
    print(f"  Plotting snapshot {n}: t = {t:.2f}")

    D = pp.pload(n, w_dir=wdir, datatype='vtk')

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f'Hall MHD - Harris Current Sheet  (t = {t:.1f})', fontsize=16)

    cs2 = 0.5
    press = cs2 * D.rho

    fields = [
        (D.rho, r'$\rho$', 'plasma'),
        (press, r'$P = c_s^2 \rho$', 'inferno'),
        (D.vx1, r'$v_x$', 'RdBu_r'),
        (D.vx2, r'$v_y$', 'RdBu_r'),
        (D.Bx1, r'$B_x$', 'RdBu_r'),
        (D.Bx2, r'$B_y$', 'RdBu_r'),
    ]

    ax_flat = axes.ravel()
    for i, (data, label, cmap) in enumerate(fields):
        ax = ax_flat[i]
        if np.abs(data).max() < 1e-14:
            continue
        if data.min() < 0 and data.max() > 0:
            vmax = max(abs(data.min()), abs(data.max()))
            norm = SymLogNorm(linthresh=vmax/100, vmin=-vmax, vmax=vmax)
        else:
            norm = None
        im = ax.pcolormesh(D.x1, D.x2, data.T, shading='auto', norm=norm, cmap=cmap)
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(label)
        ax.set_aspect('equal')

    ax_mag = ax_flat[5]
    Bmag = np.sqrt(D.Bx1**2 + D.Bx2**2)
    im = ax_mag.pcolormesh(D.x1, D.x2, Bmag.T, shading='auto', cmap='plasma')
    plt.colorbar(im, ax=ax_mag, shrink=0.8)
    ax_mag.set_xlabel('x')
    ax_mag.set_ylabel('y')
    ax_mag.set_title(r'$|\mathbf{B}|$')
    ax_mag.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(f'hall_cs_{n:04d}_t{t:.1f}.png', dpi=120)
    plt.close()

print("All plots saved!")
