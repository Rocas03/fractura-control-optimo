# Paquete `Fractura`

Todo el cómputo del proyecto. Cada número y cada figura del póster y del informe sale de aquí: los scripts
escriben los json de `../final/resultados/` y `../avance/resultados/`, y de esos json salen las cifras del texto.

## Qué hay en `src/`

Modelo y herramientas comunes:

- `modelo.jl`: los parámetros del caso base (`Par`) y el sistema adimensional en (ℓ, P).
- `integrador.jl`: Vern9 con tolerancias estrictas y un RK4 de paso fijo para el método directo y su adjunto.
- `raices.jl`: búsqueda de raíces por regula falsi (variante de Illinois).

Método directo, Pontryagin y HJB:

- `directo.jl`: control constante a trozos, costo y gradiente por adjunto, optimización con Ipopt.
- `disparo.jl`: disparo sobre el sistema extremal de Pontryagin (dos incógnitas).
- `comparaciones.jl`: inyecciones constantes y rampas contra el extremal.
- `hjb.jl`: esquema semi-Lagrangiano para la ecuación de HJB.

Un archivo por objetivo del proyecto:

- `objetivo1.jl`: dinámica sin control y con caudal constante (C1 del póster).
- `objetivo2.jl`: linealización, gramianos, observabilidad midiendo solo la presión y LQ de seguimiento (C2).
- `objetivo3.jl`: Pontryagin contra el método directo y contra controles simples (C3).
- `objetivo4.jl`: HJB, retroalimentación y lazo cerrado (C4).
- `acuerdo.jl`: el directo, el disparo y HJB en una misma tabla de convergencia (C4).
- `at2.jl`, `at2_adjunto.jl`, `objetivo5.jl`: el modelo de daño AT2 en 1D, su adjunto discreto y el control del
  programa de carga (C5).
- `sensibilidad.jl`, `fisher.jl`, `diseno.jl`: sensibilidades de la presión, matriz de Fisher y frontera entre
  costo de control e información (bloque «Diseño de Fisher»).

Verificación y salidas:

- `semilla.jl`: comprobaciones numéricas, en BigFloat, de las identidades del modelo.
- `regresion.jl`: compara cada resultado nuevo con el json que ya está en disco antes de sobrescribirlo.
- `cifras.jl`, `figuras.jl`, `datos_beamer.jl`: cifras, figuras y curvas del reporte de avance.
- `cifras_informe.jl`: cifras del informe final, una macro de LaTeX por valor de los json.

## Scripts

Se corren desde la raíz del repositorio, con Julia 1.11.4 (ver el README principal):

- `test/runtests.jl`: las pruebas.
- `regenerar.jl`: resultados de la entrega de avance (objetivos 1, 3 y 4) y `../semilla/`; además, cifras,
  figuras y curvas del reporte y del video de avance, en `../avance/reporte/` y `../avance/video/`.
- `regenerar_final.jl`: resultados y figuras de la entrega final (objetivos 2 y 5, acuerdo de métodos y
  diseño de Fisher); incluye `cruce_atribucion.jl`.
- `bin/jl`: envoltorio que llama al binario de Julia 1.11.4 en la ruta donde lo instala `juliaup` en un Mac
  con Apple Silicon; para otra ruta, se pasa en la variable `JULIA_BIN`.

Los comentarios del código citan a veces `docs/` o el reporte de avance: son notas de trabajo del proyecto que
no forman parte de este repositorio.
