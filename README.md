# 🌌 Astrofísica Computacional — 2026-I

![SUNdead](funyanimation.gif)

---

### 👨‍🏫 Profesor: José Iván Campos
### 🎓 Estudiante: Carlos Alberto Martínez Sibaja
### 📚 Curso: Astrofísica Computacional · Semestre 2026-I

---

## Contenido del repositorio

```
astrofisica-computacional/
│
├── tareas/                        ← Tareas del curso
│   ├── tareas_resueltas.md        ← Índice unificado con resúmenes
│   └── notebooks/                 ← Jupyter notebooks de cada tarea (01–15)
│
├── proyecto/                      ← Proyecto 1: Viento Solar de Parker
│   ├── README.md                  ← Documentación completa del proyecto
│   ├── constantes.py              ← Constantes físicas (SI)
│   ├── analitico.py               ← Solución semi-analítica (búsqueda de raíces)
│   ├── solucionadores_numericos.py ← Integradores Euler y RK4 (sin scipy ODE)
│   ├── generar_notebook.py        ← Constructor programático del notebook
│   └── viento_solar_parker.ipynb  ← Reporte principal
│
├── proyecto2/                     ← Proyecto 2: Reconexión Magnética Hall MHD
│   ├── README.md                  ← Documentación completa del proyecto
│   ├── Current_Sheet/             ← Lámina de corriente de Harris con PLUTO
│   │   ├── init.c                 ← Condición inicial de Harris
│   │   ├── definitions_01.h       ← Configuración Hall MHD
│   │   ├── pluto_01.ini           ← Parámetros de la corrida
│   │   ├── analysis/              ← Postproceso y figuras
│   │   ├── python_reproduction/   ← Solver Hall-MHD 2.5D en Python
│   │   └── report/report.md       ← Reporte técnico
│   └── Whistler_Waves/            ↑ Configuraciones para ondas whistler
│
└── README.md                      ← Estás aquí
```

---

## 📝 Tareas

15 tareas resueltas que cubren los algoritmos fundamentales del curso:

| Bloque | Temas |
|--------|-------|
| **Diferenciación** | Diferencias finitas, convergencia, cancelación catastrófica |
| **Integración** | Riemann, trapezoide, Simpson 1/3 y 3/8 |
| **Búsqueda de raíces** | Bisección, Newton-Raphson, secante |
| **Optimización** | Sección áurea, máximos de funciones unimodales |
| **EDOs** | Euler, RK2 (Ralston), RK4 (Runge y Kutta) |
| **Aplicaciones** | Galaxia espiral, Betelgeuse, cohete, decaimiento C-14, viento solar |

👉 Ver el índice completo: [`tareas/tareas_resueltas.md`](tareas/tareas_resueltas.md)

---

## 🚀 Proyecto: Viento Solar de Parker

Simulación numérica del viento solar transónico resolviendo la EDO de Parker
con integradores codificados a mano (Euler + RK4), validados contra la solución
semi-analítica por búsqueda de raíces (método de Brent).

**Resultados clave** ($T = 10^6$ K, $\mu = 0.5$):

| Cantidad | Valor |
|----------|-------|
| Velocidad crítica $v_c$ | ~91 km/s |
| Radio crítico $r_c$ | ~5.8 R☉ |
| Velocidad en 1 UA | ~400 km/s |
| Error relativo RK4 en 50 R☉ | ~10⁻⁶ |

👉 Documentación completa: [`proyecto/README.md`](proyecto/README.md)

---

## 🚀 Proyecto 2: Reconexión Magnética Hall MHD

Simulación de reconexión magnética en una lámina de corriente de Harris usando
Hall MHD con el código PLUTO. Incluye un solver educativo Hall-MHD 2.5D en Python
con hiperdisipación de 4to orden y diagnósticos cuantitativos (flujo reconectado,
corriente máxima, tasa de reconexión).

**Resultados clave** (malla 256×128, $t = 60$):

| Cantidad | Valor |
|----------|-------|
| Flujo reconectado | 4.55 |
| Corriente máxima $\max\|J_z\|$ | 3.08 (cerca de $t = 25$) |
| Tasa de reconexión máxima | 0.18 |

👉 Documentación completa: [`proyecto2/README.md`](proyecto2/README.md)

---

*Repositorio activo · Última actualización: junio 2026*
