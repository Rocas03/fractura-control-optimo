# Control óptimo de la tasa de inyección en un modelo reducido de fractura hidráulica

Proyecto de **MA4703-1 Control Óptimo: Teoría y Laboratorio** (FCFM, Universidad de Chile, Primavera 2026).

**Autores:** Cristóbal Castellanos y Esteban Ibarra.

El modelo reduce la fractura a dos estados, el largo $\ell$ y la presión $P$, y toma como control el caudal
inyectado $U\in[0,q_{\max}]$. Sobre él se estudian la dinámica sin control, el control LQ del sistema
linealizado, el principio de Pontryagin, la ecuación de Hamilton-Jacobi-Bellman y su comparación con un
método directo y uno de disparo, un problema de control sobre un modelo phase-field AT2 y un diseño de
experimento basado en la matriz de Fisher.

## Qué hay

| Ruta | Qué es |
|---|---|
| `julia/` | Todo el cómputo, como paquete `Fractura`: `src/` (una sección por archivo), `test/`, `regenerar.jl`, `regenerar_final.jl`, y `Project.toml`/`Manifest.toml` con las versiones exactas |
| `notebook/` | `fractura.ipynb`, notebook autocontenido, y `derivar.py`, que lo arma desde `julia/src/` |
| `final/resultados/`, `final/figuras/` | Resultados (json) y figuras (pgfplots) del informe final, escritos por `julia/regenerar_final.jl` |
| `avance/resultados/` | Resultados de la etapa de avance; las pruebas y `cifras_informe.jl` los leen como referencia |
| `semilla/` | Comprobaciones de lectura del modelo, escritas por `julia/regenerar.jl` |

El informe final y el póster se agregan después de su entrega.

## Cómo correr las pruebas

El `Manifest.toml` está fijado a **Julia 1.11.4**. Con esa versión:

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia julia/test/runtests.jl
```

`julia/bin/jl` es un envoltorio que apunta al binario 1.11.4 instalado por `juliaup`; si tu ruta es otra,
pásala en la variable `JULIA_BIN`. `final/regenerar.sh` regenera los json y las figuras de `final/` y además
necesita `latexmk`, `xelatex` y `pdfinfo`.

## Licencia

MIT, ver [`LICENSE`](LICENSE).
