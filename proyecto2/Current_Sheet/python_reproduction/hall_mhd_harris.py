"""
hall_mhd_harris.py
Reproduccion Python - Harris Current Sheet.

Este script implementa el pipeline Python auxiliar del problema:
  1. Setup de la condicion inicial (Harris sheet)
  2. Calculo de cantidades derivadas (J, flujo reconectado, etc.)
  3. Visualizacion (comparacion con PLUTO)
  4. Analisis de errores

Las ecuaciones MHD y Hall se resuelven con PLUTO. Este script
reproduce el setup y procesa los datos de PLUTO para analisis y visualizacion.
"""
from pathlib import Path
import os
import sys
import warnings

MPLCONFIG = Path('/tmp') / 'matplotlib-cache'
MPLCONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(MPLCONFIG))

import numpy as np
from numpy import pi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / 'output'
PLUTO_RUN = ROOT / 'pluto_sim'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. SETUP: Condicion inicial de Harris sheet
# ============================================================
Lx, Ly = 25.6, 12.8
nx, ny = 256, 128
dx, dy = Lx/nx, Ly/ny

xc = np.linspace(-Lx/2+dx/2, Lx/2-dx/2, nx)
yc = np.linspace(-Ly/2+dy/2, Ly/2-dy/2, ny)
X, Y = np.meshgrid(xc, yc, indexing='ij')

cs2 = 0.5; cs = np.sqrt(cs2)
B0 = 1.0; l = 0.5; Psi0 = 0.02

rho0 = 0.2 + 1.0/np.cosh(Y/l)**2
vx0 = np.zeros((nx, ny))
vy0 = np.zeros((nx, ny))
Bx0 = B0 * np.tanh(Y/l)
By0 = np.zeros((nx, ny))

kx = pi/Lx; ky = pi/Ly
Bx0 += -Psi0*ky*np.sin(ky*Y)*np.cos(2*kx*X)
By0 += Psi0*2*kx*np.sin(2*kx*X)*np.cos(ky*Y)

# ============================================================
# 2. FUNCIONES DE ANALISIS
# ============================================================
def compute_current(Bx, By, dx, dy):
    """Corriente Jz = dBy/dx - dBx/dy."""
    Jz = np.zeros_like(Bx)
    Jz[1:-1,:] = (By[2:,:] - By[:-2,:])/(2*dx)
    Jz[:,1:-1] -= (Bx[:,2:] - Bx[:,:-2])/(2*dy)
    return Jz

def reconnected_flux(By, xc, yc):
    """Flujo reconectado usado en el reporte: int_0^(Lx/2) |By(x,0)| dx."""
    j0 = int(np.argmin(np.abs(yc)))
    xmax = 0.5 * (xc.max() - xc.min())
    mask = (xc >= 0.0) & (xc <= xmax)
    return np.trapezoid(np.abs(By[mask, j0]), xc[mask])

def compute_divB(Bx, By, dx, dy):
    """Divergencia de B."""
    divB = np.zeros_like(Bx)
    divB[1:-1,:] = (Bx[2:,:] - Bx[:-2,:])/(2*dx)
    divB[:,1:-1] += (By[:,2:] - By[:,:-2])/(2*dy)
    return divB

# ============================================================
# 3. CARGAR DATOS DE PLUTO
# ============================================================
pluto_py_candidates = []
if os.environ.get('PLUTO_DIR'):
    pluto_py_candidates.append(Path(os.environ['PLUTO_DIR']) / 'Tools' / 'pyPLUTO')
pluto_py_candidates.append(ROOT.parent / 'PLUTO' / 'Tools' / 'pyPLUTO')
for candidate in pluto_py_candidates:
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))
        break
import pyPLUTO.pload as pp

with open(PLUTO_RUN / 'vtk.out') as f:
    lines = f.readlines()

snapshots = []
for line in lines:
    parts = line.split()
    snapshots.append({'n': int(parts[0]), 'time': float(parts[1])})

print(f"PLUTO: {len(snapshots)} snapshots cargados")

pluto = []
for s in snapshots:
    D = pp.pload(s['n'], w_dir=str(PLUTO_RUN) + '/', datatype='vtk')
    pluto.append({
        't': s['time'], 'n': s['n'],
        'rho': D.rho, 'vx': D.vx1, 'vy': D.vx2,
        'Bx': D.Bx1, 'By': D.Bx2,
        'x1': D.x1, 'x2': D.x2
    })
    print(f"  t={s['time']:.1f}: OK")

# ============================================================
# 4. ANALISIS: Flujo reconectado, corriente, etc.
# ============================================================
times = [d['t'] for d in pluto]
flux_rec = [reconnected_flux(d['By'], d['x1'], d['x2']) for d in pluto]
Jz_max = [np.max(np.abs(compute_current(d['Bx'], d['By'], dx, dy))) for d in pluto]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
ax1.plot(times, flux_rec, 'o-', color='C0')
ax1.set_xlabel('Time'); ax1.set_ylabel('Reconnected Flux')
ax1.set_title('Flujo Reconectado vs Tiempo (PLUTO Hall MHD)')
ax1.grid(True)

ax2.plot(times, Jz_max, 's-', color='C3')
ax2.set_xlabel('Time'); ax2.set_ylabel('Max |Jz|')
ax2.set_title('Corriente Maxima vs Tiempo')
ax2.grid(True)
plt.tight_layout()
plt.savefig(OUT_DIR / 'analysis_timeseries.png', dpi=120)
plt.close()
print("Analisis temporal guardado.")

# ============================================================
# 5. COMPARACION: Estado final PLUTO vs inicial
# ============================================================
Df = pluto[-1]
D0 = pluto[0]

fig, axes = plt.subplots(4, 3, figsize=(18, 20))
fig.suptitle('Comparacion PLUTO: Inicial (t=0) vs Final (t=60)', fontsize=16)

variables = [
    (D0['rho'], Df['rho'], r'$\rho$', 'plasma', False),
    (cs2*D0['rho'], cs2*Df['rho'], r'$P$', 'inferno', False),
    (D0['vx'], Df['vx'], r'$v_x$', 'RdBu_r', True),
    (D0['vy'], Df['vy'], r'$v_y$', 'RdBu_r', True),
]

for row, (v0, vf, label, cmap, sym) in enumerate(variables):
    for col, (data, title) in enumerate([(v0, 't=0'), (vf, f't={Df["t"]:.1f}')]):
        ax = axes[row, col]
        nrm = None
        if sym and np.max(np.abs(data)) > 1e-14:
            vm = max(abs(data.min()),abs(data.max()))
            nrm = SymLogNorm(linthresh=vm/100, vmin=-vm, vmax=vm)
        im = ax.pcolormesh(Df['x1'], Df['x2'], data.T, shading='auto', norm=nrm, cmap=cmap)
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_xlabel('x'); ax.set_ylabel('y')
        ax.set_title(f'{label} - {title}'); ax.set_aspect('equal')

# Fila de B
Bmag0 = np.sqrt(D0['Bx']**2 + D0['By']**2)
Bmagf = np.sqrt(Df['Bx']**2 + Df['By']**2)
for col, (data, title) in enumerate([(Bmag0, 't=0'), (Bmagf, f't={Df["t"]:.1f}')]):
    ax = axes[3, col]
    im = ax.pcolormesh(Df['x1'], Df['x2'], data.T, shading='auto', cmap='plasma')
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'$|\\mathbf{{B}}|$ - {title}'); ax.set_aspect('equal')

axes[3, 2].axis('off')
plt.tight_layout()
plt.savefig(OUT_DIR / 'comparison_initial_final.png', dpi=120)
plt.close()
print("Comparacion inicial/final guardada.")

# ============================================================
# 6. ERRORES y divergencia
# ============================================================
divB = [compute_divB(d['Bx'], d['By'], dx, dy) for d in pluto]
divB_L2 = [np.sqrt(np.mean(d**2)) for d in divB]

fig, ax = plt.subplots(figsize=(8, 5))
ax.semilogy(times, divB_L2, '^-', color='purple')
ax.set_xlabel('Time'); ax.set_ylabel('L2 norm of div(B)')
ax.set_title('Evolucion de div(B) (error de divergence cleaning)')
ax.grid(True)
plt.tight_layout()
plt.savefig(OUT_DIR / 'divB_evolution.png', dpi=120)
plt.close()
print("div(B) evolution guardada.")

# ============================================================
# 7. PERFILES 1D (corte en x=0)
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('Perfiles 1D en x=0: PLUTO Hall MHD', fontsize=16)
ix = nx // 2

for i, (var, label, data_list) in enumerate([
    (lambda d: d['rho'], r'$\rho$', pluto),
    (lambda d: cs2*d['rho'], r'$P$', pluto),
    (lambda d: d['vx'], r'$v_x$', pluto),
    (lambda d: d['vy'], r'$v_y$', pluto),
    (lambda d: d['Bx'], r'$B_x$', pluto),
    (lambda d: d['By'], r'$B_y$', pluto),
]):
    ax = axes.flat[i]
    for d in pluto[::3]:
        ax.plot(d['x2'], var(d)[ix,:], label=f't={d["t"]:.0f}', alpha=0.7)
    ax.set_xlabel('y'); ax.set_title(label)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_DIR / 'profiles_x0.png', dpi=120)
plt.close()
print("Perfiles 1D guardados.")

# ============================================================
# 8. RESUMEN
# ============================================================
print("\n" + "="*50)
print("RESUMEN: Harris Current Sheet con Hall MHD")
print("="*50)
print(f"Dominio: [{xc.min():.1f}, {xc.max():.1f}] x [{yc.min():.1f}, {yc.max():.1f}]")
print(f"Grilla: {nx} x {ny}")
print(f"Ancho lamina: l = {l}")
print(f"Perturbacion: Psi0 = {Psi0}")
print(f"Tiempo final: t = {Df['t']:.1f}")
print(f"Pasos PLUTO: {len(snapshots)} snapshots")
print(f"Flujo reconectado final: {flux_rec[-1]:.4f}")
print(f"|B| max final: {np.max(Bmagf):.4f}")
print(f"Archivos de salida en {OUT_DIR}")
