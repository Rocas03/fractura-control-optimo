# Control óptimo de la tasa de inyección en un modelo reducido de fractura hidráulica

<h2 align="center"><a href="https://claude.ai/artifact/1Df9xkG6aArXdK6PZSyjto">▶ Abre el póster vivo</a></h2>
<p align="center"><b><a href="https://claude.ai/artifact/1Df9xkG6aArXdK6PZSyjto">Toca aquí</a></b>: cambia el problema y mira cómo cambian todos los gráficos del póster.</p>

---

Proyecto de MA4703-1 Control Óptimo: Teoría y Laboratorio (FCFM, Universidad de Chile, primavera 2026),
profesor Héctor Ramírez C.

Cristóbal Castellanos y Esteban Ibarra.

## Si llegaste por el QR del póster

- [Póster](entrega/poster-MA4703-Castellanos-Ibarra.pdf), en la versión A0 entregada.
- [Informe final](entrega/informe-final.pdf), 15 páginas: el modelo, las demostraciones, todas las cifras y lo
  que queda abierto.
- [Notebook del avance](notebook/fractura.ipynb), que GitHub muestra con sus salidas y figuras.
- [Código](julia/): el paquete de Julia que calcula cada número del póster y del informe.

## El problema

La fractura se reduce a dos estados adimensionales, el largo ℓ y la presión P, y el control es el caudal
inyectado U, con 0 ≤ U ≤ q<sub>max</sub>. Es el modelo KGD en la forma de Geertsma y de Klerk [5], en el
régimen dominado por la tenacidad [4], con una ley de propagación subcrítica de tipo potencia [1] y una
filtración hacia la roca de coeficiente η.

Se busca un programa de inyección que, en un horizonte T, lleve la fractura a un largo objetivo ℓ<sub>d</sub>
manteniendo baja la presión y gastando poco fluido. El costo J suma el cuadrado del error de largo en T y la integral de los
cuadrados de P y de U, con pesos w<sub>T</sub>, w<sub>P</sub> y w<sub>q</sub>. Los parámetros son sintéticos y
adimensionales; no hay calibración con datos reales.

## Qué se obtuvo

Los títulos siguen los bloques del póster.

**Dinámica (C1).** Con filtración positiva y sin inyectar, cada trayectoria termina detenida y despresurizada
en un largo límite ℓ<sub>∞</sub>. Con caudal constante positivo nunca hay equilibrio.

**Control LQ (C2).** El sistema linealizado en torno a la trayectoria nominal es controlable y, midiendo solo
la presión, observable. Si la presión inicial o η son un 20 % mayores que los nominales, el regulador LQ de
seguimiento reduce la desviación terminal en un factor 15,2 y 3,3, respectivamente, frente al lazo abierto.

**Pontryagin (C3).** El control óptimo es U\* = proy<sub>[0, q<sub>max</sub>]</sub>(ζ<sub>P</sub> / (w<sub>q</sub> D)),
donde ζ<sub>P</sub> es el precio sombra de la presión y D = c + ℓ² el almacenamiento. Por transversalidad
ζ<sub>P</sub>(T) = 0, así que el caudal óptimo se anula al final. En el caso base la cota q<sub>max</sub> no se
activa y el extremal concentra el caudal al inicio: la mejor rampa decreciente queda a 0,11 % de su costo, la
mejor inyección constante a 6,7 %, y una constante con el mismo volumen cuesta 22 % más.

**HJB y acuerdo de métodos (C4).** La ecuación de HJB se resuelve con un esquema semi-Lagrangiano en una malla
de (ℓ, P, τ), y de ella sale una retroalimentación. El método directo (Ipopt, N = 200), el disparo sobre el
sistema de Pontryagin y el lazo cerrado de HJB dan el mismo costo, con diferencia relativa máxima 1,4 · 10⁻⁴.

**Phase-field AT2 (C5).** El problema se traslada a un modelo de daño por gradiente AT2 en una barra 1D, con
un defecto de G<sub>c</sub> y daño irreversible. Ahí se controla el programa de carga, no el caudal. El
descenso proyectado, con el gradiente calculado por un adjunto discreto, queda a 0,53 % del óptimo cerrado.

**Diseño de Fisher.** Con 40 medidas de presión con ruido σ = 0,001, se agrega al costo un término
−γ log det M, con M la matriz de información de Fisher. En el codo de la frontera (γ = 0,03), la cota de
Cramér-Rao para η baja 2,6 veces a cambio de 6 % más de costo. No se afirma identificabilidad.

**Lo que queda abierto.** Los tres métodos coinciden, pero no está probado que el extremal calculado sea el
minimizador. El modelo es reducido: entre él y el modelo phase-field hay una analogía, no un teorema.

## Dónde está cada cosa

- `entrega/`: el póster y el informe final, en PDF.
- `julia/`: todo el cómputo, como paquete de Julia `Fractura`. Su [README](julia/README.md) dice qué
  archivo hace qué.
- `notebook/`: el notebook de la entrega de avance (26 de octubre), autocontenido; calcula los objetivos 1, 3
  y 4 y corre sus pruebas. `derivar.py` lo arma desde `julia/src/`, pero necesita el reporte de avance, que no
  está en este repositorio.
- `final/resultados/` y `final/figuras/`: los json y las figuras (pgfplots) del informe final.
- `avance/resultados/`: los json de la entrega de avance, que las pruebas usan como referencia.
- `semilla/`: comprobaciones numéricas de la lectura del modelo.

## Cómo reproducir los resultados

El `Manifest.toml` está fijado a Julia 1.11.4. Con [juliaup](https://github.com/JuliaLang/juliaup), desde la
raíz del repositorio:

```sh
juliaup add 1.11.4
julia +1.11.4 --project=julia -e 'using Pkg; Pkg.instantiate()'
julia +1.11.4 --project=julia julia/test/runtests.jl
```

La primera vez baja y precompila los paquetes. Con los paquetes ya instalados, las 720 pruebas pasan en menos
de un minuto (medido el 24 de noviembre de 2026).

Para regenerar los resultados:

```sh
julia +1.11.4 --project=julia julia/regenerar_final.jl
julia +1.11.4 --project=julia julia/regenerar.jl
```

El primero escribe los json y las figuras de `final/` y tarda unos cuatro minutos; el segundo, los json de la
entrega de avance y de `semilla/`. En una copia limpia del repositorio, los dos reproducen byte a byte los
archivos versionados (medido el 24 de noviembre de 2026). También escriben las cifras del informe y los archivos del
reporte y del video de avance, en `final/informe/`, `avance/reporte/` y `avance/video/`, carpetas que no están
en el repositorio.

## Referencias

Son las del informe final. El póster cita [1], [2], [5], [7] y [8].

1. B. K. Atkinson. Subcritical crack growth in geological materials. *Journal of Geophysical Research: Solid
   Earth* 89(B6):4077–4114, 1984. [doi:10.1029/JB089iB06p04077](https://doi.org/10.1029/JB089iB06p04077)
2. B. Bourdin, G. A. Francfort, J.-J. Marigo. *The Variational Approach to Fracture*. Springer, 2008.
   [doi:10.1007/978-1-4020-6395-4](https://doi.org/10.1007/978-1-4020-6395-4)
3. B. Bourdin, A. Peirce. *Modeling seismicity induced by hydraulic fracturing in the El Teniente copper
   mine*. White paper no publicado, 20 de diciembre de 2018.
4. E. Detournay. Propagation regimes of fluid-driven fractures in impermeable rocks. *International Journal
   of Geomechanics* 4(1):35–45, 2004.
   [doi:10.1061/(ASCE)1532-3641(2004)4:1(35)](https://doi.org/10.1061/(ASCE)1532-3641(2004)4:1(35))
5. J. Geertsma, F. de Klerk. A rapid method of predicting width and extent of hydraulically induced
   fractures. *Journal of Petroleum Technology* 21(12):1571–1581, 1969.
   [doi:10.2118/2458-PA](https://doi.org/10.2118/2458-PA)
6. Q. Gu, K. A. Hoo. Model-based closed-loop control of the hydraulic fracturing process. *Industrial &
   Engineering Chemistry Research* 54(5):1585–1594, 2015.
   [doi:10.1021/ie5024782](https://doi.org/10.1021/ie5024782)
7. I. Neitzel, T. Wick, W. Wollner. An optimal control problem governed by a regularized phase-field fracture
   propagation model. *SIAM Journal on Control and Optimization* 55(4):2271–2288, 2017.
   [doi:10.1137/16M1062375](https://doi.org/10.1137/16M1062375)
8. E. Trélat. *Contrôle optimal : théorie et applications*. Versión electrónica, 2013 (1.ª edición: Vuibert,
   Collection Mathématiques Concrètes, 2005).
   [ljll.fr/trelat/fichiers/livreopt.pdf](https://www.ljll.fr/trelat/fichiers/livreopt.pdf)
9. Y. Navidtehrani, C. Betegón, J. Vallejos, E. Martínez-Pañeda. A phase field model for hydraulic fracture:
   Drucker–Prager driving force and a hybrid coupling strategy. *Computer Methods in Applied Mechanics and
   Engineering* 444:118155, 2025. [doi:10.1016/j.cma.2025.118155](https://doi.org/10.1016/j.cma.2025.118155)

## Licencia

MIT; ver [`LICENSE`](LICENSE).
