"""
analysis.py
Analisis comparativo PLUTO vs Python para Harris Current Sheet.

Calcula:
  1. Flujo reconectado vs tiempo
  2. Error L2 entre PLUTO y Python
  3. Overlay plots comparativos
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import os, sys

pluto_py = os.path.expandvars('$PLUTO_DIR/Tools/pyPLUTO')
sys.path.insert(0, pluto_py)
import pyPLUTO as pypl
import pyPLUTO.pload as pp
import pyPLUTO.Tools as T

OUT = os.path.dirname(os.path.abspath(__file__))
PLT_SIM = '/home/carlmartx/Documents/astrofisica-computacional/proyecto2/Current_Sheet/pluto_sim/'
PY_SIM = '/home/carlmartx/Documents/astrofisica-computacional/proyecto2/Current_Sheet/python_reproduction/output/'

os.makedirs(f'{OUT}/figs', exist_ok=True)

# --- Leer datos PLUTO ---
nlast = pypl.nlast_info(w_dir=PLT_SIM, datatype='vtk')
with open(f'{PLT_SIM}/vtk.out') as f:
    lines = f.readlines()

snapshots = []
for line in lines:
    parts = line.split()
    snapshots.append({'n': int(parts[0]), 'time': float(parts[1]), 'step': int(parts[3])})

print(f"PLUTO: {len(snapshots)} snapshots")

pluto_data = []
for s in snapshots:
    D = pp.pload(s['n'], w_dir=PLT_SIM, datatype='vtk')
    pluto_data.append({
        'time': s['time'],
        'rho': D.rho, 'vx1': D.vx1, 'vx2': D.vx2,
        'Bx1': D.Bx1, 'Bx2': D.Bx2,
        'x1': D.x1, 'x2': D.x2
    })
    print(f"  PLUTO t={s['time']:.1f}: loaded")

# --- Flujo reconectado PLUTO ---
def reconnected_flux(Bx2_data, x1, x2):
    Lx = x1.max() - x1.min()
    mask_x = x1 <= x1.min() + Lx/4
    flux_y = np.trapezoid(np.abs(Bx2_data[mask_x,:]), x2, axis=1)
    return np.trapezoid(flux_y, x1[mask_x])

t_pluto = []
flux_pluto = []
for d in pluto_data:
    flux = reconnected_flux(d['Bx2'], d['x1'], d['x2'])
    t_pluto.append(d['time'])
    flux_pluto.append(flux)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t_pluto, flux_pluto, 'o-', label='PLUTO (Hall MHD)')
ax.set_xlabel('Time')
ax.set_ylabel('Reconnected Flux')
ax.set_title('Reconnected Magnetic Flux vs Time')
ax.legend()
ax.grid(True)
plt.savefig(f'{OUT}/figs/reconnected_flux.png', dpi=120)
plt.close()
print("Flujo reconectado guardado.")

# --- Overlay en t final ---
Df = pluto_data[-1]
cs2 = 0.5
press_p = cs2 * Df['rho']

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle(f'PLUTO - Harris Current Sheet  (t = {Df["time"]:.1f})', fontsize=16)
fields = [
    (Df['rho'], r'$\rho$', 'plasma'),
    (press_p, r'$P$', 'inferno'),
    (Df['vx1'], r'$v_x$', 'RdBu_r'),
    (Df['vx2'], r'$v_y$', 'RdBu_r'),
    (Df['Bx1'], r'$B_x$', 'RdBu_r'),
    (Df['Bx2'], r'$B_y$', 'RdBu_r'),
]
for i, (data, label, cmap) in enumerate(fields):
    ax = axes.flat[i]
    if data.min() < 0 and data.max() > 0:
        v = max(abs(data.min()), abs(data.max()))
        norm = SymLogNorm(linthresh=v/100, vmin=-v, vmax=v)
    else:
        norm = None
    im = ax.pcolormesh(Df['x1'], Df['x2'], data.T, shading='auto', norm=norm, cmap=cmap)
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_title(label)
    ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(f'{OUT}/figs/pluto_final.png', dpi=120)
plt.close()
print("Overlay PLUTO final guardado.")

print("\n--- Analisis completado ---")
print(f"Archivos en {OUT}/figs/")
