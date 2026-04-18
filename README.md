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
│   └── notebooks/                 ← Jupyter notebooks de cada tarea (01–12)
│
├── proyecto/                      ← Proyecto final: Viento Solar de Parker
│   ├── README.md                  ← Documentación completa del proyecto
│   ├── constantes.py              ← Constantes físicas (SI)
│   ├── analitico.py               ← Solución semi-analítica (búsqueda de raíces)
│   ├── solucionadores_numericos.py ← Integradores Euler y RK4 (sin scipy ODE)
│   ├── generar_notebook.py        ← Constructor programático del notebook
│   └── viento_solar_parker.ipynb  ← Reporte principal
│
└── README.md                      ← Estás aquí
```

---

## 📝 Tareas

12 tareas resueltas que cubren los algoritmos fundamentales del curso:

| Bloque | Temas |
|--------|-------|
| **Diferenciación** | Diferencias finitas, convergencia, cancelación catastrófica |
| **Integración** | Riemann, trapezoide, Simpson 1/3 y 3/8 |
| **Búsqueda de raíces** | Bisección, Newton-Raphson, secante |
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

*Repositorio activo · Última actualización: abril 2026*
