# Proyecto 2 - Astrofisica Computacional

Este repositorio contiene la simulacion y el analisis de reconexion magnetica en una lamina de corriente de Harris usando Hall MHD con PLUTO. El caso principal esta en `Current_Sheet/`.

## Estructura

```text
.
├── Current_Sheet/
│   ├── init.c                         # Condicion inicial de Harris para PLUTO
│   ├── definitions_01.h               # Configuracion Hall MHD principal
│   ├── definitions_02.h               # Configuracion alternativa
│   ├── pluto_01.ini                   # Parametros de la corrida principal
│   ├── pluto_02.ini                   # Parametros alternativos
│   ├── pluto_sim/                     # Metadatos y salidas de la corrida PLUTO
│   ├── analysis/                      # Postproceso de snapshots PLUTO
│   │   ├── analysis.py                # Diagnosticos cuantitativos y figuras finales
│   │   ├── plot_results.py            # Paneles por snapshot
│   │   └── figs/                      # Figuras y CSV usados en el reporte
│   ├── plots/                         # Paneles temporales t=0,...,60
│   ├── python_reproduction/           # Simulacion Python Hall-MHD y diagnosticos
│   │   ├── hall_mhd_harris.py         # Solver Hall-MHD 2.5D autocontenido
│   │   └── output/                    # Snapshots, GIF, figuras y CSV
│   └── report/report.md               # Reporte tecnico en Markdown
├── PLUTO/                             # Codigo PLUTO local usado como dependencia externa
└── project_2.pdf                      # Enunciado o material base del proyecto
```

## Flujo de trabajo

Desde la raiz del repositorio:

```bash
# Figuras por snapshot desde las salidas VTK de Current_Sheet/pluto_sim
python Current_Sheet/analysis/plot_results.py

# Diagnosticos, tablas CSV y figuras finales del reporte
python Current_Sheet/analysis/analysis.py

# Simulacion Python independiente en una malla reducida rapida (64x32, t=5)
cd Current_Sheet/python_reproduction
python hall_mhd_harris.py

# Simulacion Python estable documentada en el reporte (128x64, t=25)
python hall_mhd_harris.py \
  --nx 128 --ny 64 \
  --tstop 25 --output-dt 5 \
  --cfl 0.10 --cfl-hall 0.03 \
  --eta 1.0e-2 --nu 5.0e-3 \
  --eta-h 3.0e-4 --nu-h 1.0e-4

# Simulacion Python con la misma malla, tiempo y cadencia de PLUTO
python hall_mhd_harris.py --pluto-grid
```

Para repetir la corrida con PLUTO, define `PLUTO_DIR` apuntando a la instalacion local de PLUTO y usa los archivos `Current_Sheet/definitions_01.h`, `Current_Sheet/pluto_01.ini` e `Current_Sheet/init.c`. El reporte incluye el bloque completo de comandos y la interpretacion fisica de los resultados.

## Resultados principales

- La corrida Hall MHD llega hasta `t=60` en una malla `256 x 128`.
- El flujo reconectado usado como diagnostico global alcanza `4.5525`.
- La corriente maxima llega a `max |J_z| = 3.0819` cerca de `t=25`.
- La simulacion Python ahora evoluciona el mismo setup con un solver Hall-MHD 2.5D autocontenido con hiperdisipacion de 4to orden para estabilidad con baja difusion. La corrida documentada usa `128 x 64`, `tstop=25` y salidas cada `5`.
- La simulacion Python genera automaticamente snapshots `.npz`, paneles 2D, diagnosticos CSV, series temporales y un GIF animado de la evolucion en `output/evolution.gif`.

## Notas sobre el solver Python

El script `hall_mhd_harris.py` implementa un solver educativo de diferencias finitas centradas. A diferencia de PLUTO (Godunov HLL), este esquema requiere disipacion artificial para estabilidad. Se soluciono agregando **hiperdisipacion de 4to orden** ($\nu_h\nabla^4$) que disipa solo escalas de grilla sin suprimir la reconexion fisica. Ver seccion 5.2.1 del reporte para detalles.
