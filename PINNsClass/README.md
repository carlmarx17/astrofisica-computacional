# PINNs y sismología coronal — clase práctica

Material de la clase de Physics-Informed Neural Networks para Astrofísica
Computacional 2026-I.

## Idea del material

El caso de estudio es la oscilación transversal amortiguada de un lazo coronal
(modo kink). El argumento que atraviesa cuaderno y presentación:

- para el **problema directo** (ecuación y parámetros conocidos) una PINN es peor
  que RK4 por varios órdenes de magnitud, y el material lo mide en vez de
  esconderlo;
- para el **problema inverso y la asimilación de datos** (observaciones escasas,
  ruidosas y con huecos de cobertura) la física dentro de la pérdida es lo que
  convierte un ajuste en una medición.

El juicio no se hace mirando la curva, sino auditando observables físicos:
velocidad, residuo de la EDO, energía y la ley de disipación `Ė = −2βv² ≤ 0`.

## Arquitectura

Hay dos scripts fuente y un flujo en una sola dirección. **Los `.ipynb`, `.pptx`
y `.pdf` son productos: no se editan a mano.**

```
generar_cuaderno_oscilador_pinn.py
        │  (escribe)
        ▼
oscilador_armonico_pinn_fisica_solar.ipynb
        │  (al ejecutarse produce)
        ├──► figuras/*.png, figuras/*.gif
        └──► resultados.json
                    │  (lee)
                    ▼
        construir_presentacion.py  ──►  .pptx  ──►  .pdf
```

La razón de que el cuaderno se genere desde un `.py` es que un `.ipynb` ejecutado
es un archivo enorme lleno de imágenes en base64: imposible de revisar en un diff.
El script es la fuente legible; el cuaderno es su salida.

`resultados.json` es la otra pieza deliberada: el cuaderno vuelca ahí cada cifra
que la presentación necesita, y el deck las lee. **Ninguna cifra del deck está
escrita a mano**, así que no pueden desincronizarse.

### Dentro del cuaderno

Todo el contenido se apoya en ocho funciones, encadenadas de dato a diagnóstico:

| Función | Papel |
|---|---|
| `solucion_exacta(t)` | la verdad analítica, referencia de todo |
| `rk4(f, t, y0)` | integrador clásico, el rival honesto de la PINN |
| `generar_observaciones(...)` | fabrica la serie observada: pocos puntos, ruido, hueco |
| `MLP` | la red `u_θ(t)`, con la normalización temporal dentro del `forward` |
| `perdidas(model, lam_fis, obs=...)` | las tres pérdidas: datos, física, condiciones iniciales |
| `entrenar(lam_fis, obs=...)` | Adam + L-BFGS, guardando instantáneas para la animación |
| `observables(model)` | audita el modelo: `u, v, a`, residuo, energía, energía espuria |
| `problema_inverso(obs=...)` | añade `β` y `ω₀` como variables entrenables |

Tres decisiones de diseño que sostienen el resto:

1. **`lam_fis` es el único mando del experimento.** Caja negra y PINN son la
   misma red, el mismo optimizador y los mismos datos; solo cambia ese número
   (`0` frente a `30`). Cualquier diferencia es atribuible a la pérdida.
2. **`obs=(t_tensor, u_tensor)` desacopla los datos del entrenamiento.** Permite
   cambiar las observaciones —más ruido, otro hueco, otra dinámica— sin tocar
   nada más. Es el gancho sobre el que se resuelven las tareas.
3. **El residuo va adimensionalizado** (dividido por `ω₀²`). Sin eso, los tres
   términos de la pérdida tienen unidades distintas, `λ` depende de la escala
   temporal elegida y aparece la patología de gradientes.

### Dentro del constructor de la presentación

`construir_presentacion.py` es un pequeño sistema de diseño sobre `python-pptx`:
`slide()`, `caja()`, `texto()`, `titulo()`, `vinetas()`, `codigo()`, `tabla()`,
`imagen()`. Cada diapositiva son unas pocas llamadas a esos ayudantes, así que
cambiar la paleta o la tipografía es tocar las constantes de arriba del archivo.
`codigo()` calcula su propia altura y colorea palabras clave y comentarios.

## Cómo reproducirlo

```bash
python generar_cuaderno_oscilador_pinn.py
jupyter nbconvert --to notebook --execute --inplace \
    oscilador_armonico_pinn_fisica_solar.ipynb    # ~5 min en CPU
python construir_presentacion.py
libreoffice --headless --convert-to pdf Presentacion_PINNs_Fisica_Solar.pptx
```

El orden importa. Si cambias un parámetro físico en el generador y vuelves a
correr los cuatro comandos, figuras, cifras y diapositivas quedan al día solas.

## Tareas para los estudiantes

La §11 del cuaderno tiene cuatro tareas de 5–10 minutos, cada una con su croquis
ya montado y unos pocos `TODO`:

1. barrer `λ_fis` y ver qué le pasa a la energía espuria;
2. hacer crecer el hueco de cobertura hasta que la física deje de bastar;
3. correr el problema inverso desde una suposición mala y llegar a `B` y `l/a`;
4. entrenar con la física equivocada (amortiguamiento cuadrático ajustado con un
   modelo lineal) y descubrir cómo se detecta.

Se entregan la celda ejecutada y dos o tres frases por tarea. La pregunta pesa
más que el código.

## Dependencias

`numpy`, `matplotlib`, `pandas`, `torch` (CPU basta), `nbformat`, `pillow`;
`python-pptx` solo para construir la presentación.

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install numpy matplotlib pandas nbformat pillow python-pptx
```
