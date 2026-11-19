# Control óptimo de la tasa de inyección en un modelo reducido de fractura hidráulica

Proyecto de **MA4703-1 Control Óptimo: Teoría y Laboratorio** (FCFM, Universidad de Chile, Primavera 2026).

**Autores:** Cristóbal Castellanos y Esteban Ibarra.

El modelo reduce la fractura a dos estados, el largo $\ell$ y la presión $P$, y toma como control el caudal
inyectado $U\in[0,q_{\max}]$. Sobre él se estudian la dinámica sin control, el control LQ del sistema
linealizado, el principio de Pontryagin, la ecuación de Hamilton-Jacobi-Bellman y su comparación con un
método directo y uno de disparo, un problema de control sobre un modelo phase-field AT2 y un diseño de
experimento basado en la matriz de Fisher.

## Si llegaste por el QR del póster

- **Póster:** [`entrega/poster-MA4703-Castellanos-Ibarra.pdf`](entrega/poster-MA4703-Castellanos-Ibarra.pdf) (A0, la versión entregada)
- **Informe final:** [`entrega/informe-final.pdf`](entrega/informe-final.pdf), con las demostraciones, las cifras y lo que queda abierto
- **Referencias:** al final de esta página

**El resultado en una línea:** inyectar al inicio. Una rampa decreciente de caudal queda a 0,11 % del costo del
extremal de Pontryagin; la mejor inyección constante, a 6,7 %, y a igual volumen una constante cuesta 22 % más.

| Bloque del póster | Dónde está en el código (`julia/src/`) |
|---|---|
| Dinámica (C1) | `modelo.jl`, `integrador.jl`, `objetivo1.jl` |
| Control LQ (C2) | `objetivo2.jl` |
| Pontryagin (C3) | `objetivo3.jl`, `disparo.jl` (disparo sobre el sistema de PMP), `directo.jl` (método directo con Ipopt) |
| HJB y acuerdo de métodos (C4) | `hjb.jl`, `objetivo4.jl`, `acuerdo.jl` (directo, disparo y HJB: diferencia relativa máxima 1,4 · 10⁻⁴) |
| Phase-field AT2 (C5) | `at2.jl`, `at2_adjunto.jl`, `objetivo5.jl` (descenso proyectado a 0,53 % del óptimo cerrado) |
| Diseño de Fisher | `sensibilidad.jl`, `fisher.jl`, `diseno.jl` (la cota de η baja 2,6 veces pagando +6 % de costo) |

Todas las cifras del póster y del informe las escribe el código: `julia/regenerar_final.jl` produce los json de
`final/resultados/` y las figuras de `final/figuras/`.

## Qué hay

| Ruta | Qué es |
|---|---|
| `julia/` | Todo el cómputo, como paquete `Fractura`: `src/` (una sección por archivo), `test/`, `regenerar.jl`, `regenerar_final.jl`, y `Project.toml`/`Manifest.toml` con las versiones exactas |
| `notebook/` | `fractura.ipynb`, notebook autocontenido, y `derivar.py`, que lo arma desde `julia/src/` |
| `final/resultados/`, `final/figuras/` | Resultados (json) y figuras (pgfplots) del informe final, escritos por `julia/regenerar_final.jl` |
| `avance/resultados/` | Resultados de la etapa de avance; las pruebas y `cifras_informe.jl` los leen como referencia |
| `semilla/` | Comprobaciones de lectura del modelo, escritas por `julia/regenerar.jl` |
| `entrega/` | Los PDF entregados: el póster A0 y el informe final |

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
