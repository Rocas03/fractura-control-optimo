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

## Referencias

Las del informe final (`final/informe/secciones/referencias.tex`); el póster cita [1], [2], [5], [7] y [8].

1. B. K. Atkinson, Subcritical crack growth in geological materials, *Journal of Geophysical Research: Solid Earth* 89(B6):4077–4114, 1984, doi:[10.1029/JB089iB06p04077](https://doi.org/10.1029/JB089iB06p04077).
2. B. Bourdin, G. A. Francfort, J.-J. Marigo, *The Variational Approach to Fracture*, Springer, 2008, doi:[10.1007/978-1-4020-6395-4](https://doi.org/10.1007/978-1-4020-6395-4).
3. B. Bourdin, A. Peirce, *Modeling seismicity induced by hydraulic fracturing in the El Teniente copper mine*, documento de trabajo no publicado, 20 de diciembre de 2018.
4. E. Detournay, Propagation regimes of fluid-driven fractures in impermeable rocks, *International Journal of Geomechanics* 4(1):35–45, 2004, doi:[10.1061/(ASCE)1532-3641(2004)4:1(35)](https://doi.org/10.1061/(ASCE)1532-3641(2004)4:1(35)).
5. J. Geertsma, F. de Klerk, A rapid method of predicting width and extent of hydraulically induced fractures, *Journal of Petroleum Technology* 21(12):1571–1581, 1969, doi:[10.2118/2458-PA](https://doi.org/10.2118/2458-PA).
6. Q. Gu, K. A. Hoo, Model-based closed-loop control of the hydraulic fracturing process, *Industrial & Engineering Chemistry Research* 54(5):1585–1594, 2015, doi:[10.1021/ie5024782](https://doi.org/10.1021/ie5024782).
7. I. Neitzel, T. Wick, W. Wollner, An optimal control problem governed by a regularized phase-field fracture propagation model, *SIAM Journal on Control and Optimization* 55(4):2271–2288, 2017, doi:[10.1137/16M1062375](https://doi.org/10.1137/16M1062375).
8. E. Trélat, *Contrôle optimal: théorie & applications*, Vuibert, Collection Mathématiques Concrètes, 2005.
9. Y. Navidtehrani, C. Betegón, J. Vallejos, E. Martínez-Pañeda, A phase field model for hydraulic fracture: Drucker–Prager driving force and a hybrid coupling strategy, *Computer Methods in Applied Mechanics and Engineering* 444:118155, 2025, doi:[10.1016/j.cma.2025.118155](https://doi.org/10.1016/j.cma.2025.118155).

## Licencia

MIT, ver [`LICENSE`](LICENSE).
