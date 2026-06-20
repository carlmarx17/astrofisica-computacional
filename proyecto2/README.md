# Proyecto 2 - Astrofisica Computacional

Este repositorio contiene la simulacion y el analisis de reconexion magnetica en una lamina de corriente de Harris usando Hall MHD con PLUTO. El caso principal esta en `Current_Sheet/`; `Whistler_Waves/` conserva configuraciones auxiliares para ondas whistler.

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
│   ├── python_reproduction/           # Reproduccion Python del setup y diagnosticos
│   └── report/report.md               # Reporte tecnico en Markdown
├── Whistler_Waves/                    # Configuraciones PLUTO para ondas whistler
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

# Reproduccion Python de la condicion inicial y graficas auxiliares
cd Current_Sheet/python_reproduction
python hall_mhd_harris.py
```

Para repetir la corrida con PLUTO, define `PLUTO_DIR` apuntando a la instalacion local de PLUTO y usa los archivos `Current_Sheet/definitions_01.h`, `Current_Sheet/pluto_01.ini` e `Current_Sheet/init.c`. El reporte incluye el bloque completo de comandos y la interpretacion fisica de los resultados.

## Resultados principales

- La corrida Hall MHD llega hasta `t=60` en una malla `256 x 128`.
- El flujo reconectado usado como diagnostico global alcanza `4.5525`.
- La corriente maxima llega a `max |J_z| = 3.0819` cerca de `t=25`.
- La reproduccion Python de la condicion inicial coincide con PLUTO con errores relativos L2 de orden `1e-8`.

Los productos finales mas importantes estan en `Current_Sheet/analysis/figs/` y se referencian desde `Current_Sheet/report/report.md`.
