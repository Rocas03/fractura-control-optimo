# Figuras del reporte en pgfplots, escritas desde Julia: port de figura_comparacion y figura_cota_activa
# (e31f323:avance/codigo/objetivo3_pontryagin.py) y de figura_retroalimentacion (objetivo4_hjb.py), con el estilo
# de graficos.py: letra de 7 y 8 pt, colores y trazos de matplotlib, grilla al 25 % y coma decimal. Cada figura es
# un tikzpicture con las coordenadas en línea, que el reporte incluye con \input{figuras/<n>.tex}. Los tamaños de
# eje (en pt) están calibrados para que cada figura mida lo que medía en matplotlib sin cajas sobrellenas:
# 7,0 × 1,85 in la comparación; 0,40 y 0,50 del ancho de texto × 1,75 in las otras dos.
const FIGURAS = joinpath(RAIZ, "avance", "reporte", "figuras")
const TAMANOS = (comparacion = (w = 130.9, h = 106.75, sep = 40.0), cota = (w = 168.6, h = 98.8),
                 retro = (w = 219.4, h = 99.5))

# Lo que dibujan las figuras, en vectores planos (sin soluciones de EDO adentro).
function datos_figuras(f3, f4, o4)
    p = Par()
    ts = range(0.0, p.T, length = 801)
    Y = f4.rd["sol"].(ts)
    t08, U08 = curva_extremal(f3.r08, f3.p08)
    td, Ud = curva_escalonada(f3.q08, f3.p08)
    (comparacion = [(nombre, c[1], c[2], c[3], c[4]) for (nombre, c) in f3.curvas],   # (nombre, t, U, ℓ, P)
     cota = (t = t08, U = U08, td = td, Ud = Ud, N = length(f3.q08), qmax = f3.p08.qmax, T = f3.p08.T,
             conmutaciones = collect(Float64, f3.r08["conmutaciones"])),
     retro = (t = collect(ts), U = [ustar(x2(y), ζ2(y), p) for y in Y],
              niveles = [(t = lc.t, U = lc.U, nl = f["nl"], nt = f["nt"]) for (lc, f) in zip(f4.lcs, o4["niveles"])]))
end

const COLORES = raw"""\definecolor{C0}{HTML}{1F77B4}\definecolor{C1}{HTML}{FF7F0E}\definecolor{C2}{HTML}{2CA02C}%
\definecolor{rojo}{HTML}{D62728}\definecolor{gris40}{gray}{0.4}\definecolor{gris30}{gray}{0.3}%
"""
const EJE = raw"scale only axis, axis line style={line width=0.8pt}, tick style={black, line width=0.8pt}, " *
            raw"tick align=outside, major tick length=3.5pt, xtick pos=lower, ytick pos=left, grid=major, " *
            raw"grid style={black!31, opacity=0.25, line width=0.8pt}, label style={font=\footnotesize}, " *
            raw"tick label style={font=\scriptsize}, legend style={font=\scriptsize, draw=none, fill=none}, " *
            raw"legend cell align=left, /pgf/number format/use comma, /pgf/number format/1000 sep={}, " *
            raw"scaled ticks=false, enlarge x limits=0.05, enlarge y limits=0.05, xlabel={$\tau$}"
# Trazos de matplotlib: k- C0-- C1-. C2:
const TRAZO = Dict("k-" => "black, solid", "C0--" => "C0, dashed", "C1-." => "C1, dashdotted", "C2:" => "C2, dotted",
                   "C0:" => "C0, dotted", "C1--" => "C1, dashed", "C2--" => "C2, dashed",
                   "rojo-" => "rojo, solid")

function curva(io, trazo, x, y; ancho = 1.2, extra = "")
    print(io, "\\addplot[", TRAZO[trazo], ", line width=", ancho, "pt, no marks", extra, "] table[row sep=\\\\] {x y\\\\\n")
    for (a, b) in zip(x, y)
        @printf(io, "%.6g %.6g\\\\\n", a, b)
    end
    print(io, "};\n")
end

function figura_comparacion(io, d)
    s = TAMANOS.comparacion
    print(io, COLORES, "\\begin{tikzpicture}\n\\begin{groupplot}[group style={group size=3 by 1, horizontal sep=",
          s.sep, "pt}, ", EJE, ", width=", s.w, "pt, height=", s.h, "pt]\n")
    trazos = ("k-", "C0--", "C1-.", "C2:")
    for (j, rotulo) in enumerate(("caudal \$U\$", "largo \$\\ell\$", "presión \$P\$"))
        print(io, "\\nextgroupplot[ylabel={", rotulo, "}", j == 2 ? ", legend pos=north west" : "", "]\n")
        for (c, trazo) in zip(d, trazos)
            curva(io, trazo, c[2], c[2 + j])
            j == 2 && print(io, "\\addlegendentry{", c[1], "}\n")
        end
    end
    print(io, "\\end{groupplot}\n\\end{tikzpicture}%\n")
end

function figura_cota_activa(io, d)
    s = TAMANOS.cota
    print(io, COLORES, "\\begin{tikzpicture}\n\\begin{axis}[", EJE, ", width=", s.w, "pt, height=", s.h, "pt, ",
          "ylabel={caudal \$U\$}, enlarge y limits=false, ymin=-0.05, ymax=", d.qmax + 0.25, ", legend pos=south west]\n")
    curva(io, "k-", d.t, d.U)
    print(io, "\\addlegendentry{disparo}\n")
    curva(io, "C1--", d.td, d.Ud)
    print(io, "\\addlegendentry{directo, \$N=", d.N, "\$}\n")
    print(io, "\\draw[gris40, dotted, line width=0.8pt] ({rel axis cs:0,0} |- {axis cs:0,", d.qmax,
          "}) -- ({rel axis cs:1,0} |- {axis cs:0,", d.qmax, "});\n")
    print(io, "\\node[anchor=south east, xshift=-2pt, yshift=3pt, font=\\scriptsize, text=gris30] at (axis cs:",
          d.T, ",", d.qmax, ") {\$q_{\\max}\$};\n")
    for c in d.conmutaciones
        print(io, "\\draw[rojo, dashed, line width=0.8pt] ({axis cs:", c, ",0} |- {rel axis cs:0,0}) -- ({axis cs:",
              c, ",0} |- {rel axis cs:0,1});\n")
    end
    print(io, "\\end{axis}\n\\end{tikzpicture}%\n")
end

function figura_retroalimentacion(io, d)
    s = TAMANOS.retro
    print(io, COLORES, "\\begin{tikzpicture}\n\\begin{axis}[", EJE, ", width=", s.w, "pt, height=", s.h, "pt, ",
          "ylabel={caudal \$U\$}, legend pos=south west]\n")
    curva(io, "k-", d.t, d.U; ancho = 1.8)
    print(io, "\\addlegendentry{\$U^*\$ del PMP}\n")
    for (n, trazo) in zip(d.niveles, ("C0:", "C1-.", "C2--"))
        curva(io, trazo, n.t, n.U; ancho = 1.0, extra = ", const plot")
        print(io, "\\addlegendentry{HJB, malla \$", n.nl, "^2\\times", n.nt, "\$}\n")
    end
    print(io, "\\end{axis}\n\\end{tikzpicture}%\n")
end

function escribir_figuras(d, dir = FIGURAS)
    mkpath(dir)
    for (nombre, f, x) in (("comparacion", figura_comparacion, d.comparacion), ("cota-activa", figura_cota_activa, d.cota),
                           ("retroalimentacion", figura_retroalimentacion, d.retro))
        archivo = joinpath(dir, nombre * ".tex")
        open(archivo, "w") do io
            println(io, "% Escrito por julia/src (paquete Fractura). No editar a mano.")
            f(io, x)
        end
        println("figuras: $nombre.tex, $(round(Int, filesize(archivo) / 1024)) KB")
    end
end

# ---- F2 · la figura del acuerdo entre los tres métodos -------------------------------------------
# Un solo par de ejes con los tres controles, y un recuadro interior con las discrepancias punto a
# punto en escala logarítmica. El recuadro son EJES PROPIOS, no un segundo eje y sobre los mismos:
# dos escalas y en un mismo marco es el error de lectura más común en una figura y acá no haría
# falta, porque el recuadro compara |ΔU| con |ΔU| y no con U.
# El recuadro va arriba a la derecha porque ahí no hay datos: U va de 1,08 en τ = 0 a casi 0 en
# τ = T, o sea de la esquina superior izquierda a la inferior derecha; la leyenda ocupa la inferior
# izquierda, como en figura_retroalimentacion.
# Identidad por color Y por trazo (continuo / rayado / punteado): un lector que no distingue el azul
# del verde sigue separando las tres curvas.
const TAM_ACUERDO = (w = 219.4, h = 132.0, wr = 74.0, hr = 46.0)
const EJE_CHICO = raw"scale only axis, axis line style={line width=0.6pt}, tick style={black, line width=0.6pt}, " *
                  raw"tick align=outside, major tick length=2pt, xtick pos=lower, ytick pos=left, grid=major, " *
                  raw"grid style={black!31, opacity=0.25, line width=0.6pt}, label style={font=\tiny}, " *
                  raw"tick label style={font=\tiny}, /pgf/number format/use comma, /pgf/number format/1000 sep={}, " *
                  raw"scaled ticks=false, axis background/.style={fill=white}, ymode=log"

# Escalera del directo: el control es constante en cada tramo, así que se dibuja con `const plot` y
# hay que repetir el último valor en el nodo final para que el último escalón se cierre.
escalera(t, q) = (collect(t), vcat(collect(q), q[end]))

function figura_acuerdo(io, d)
    s = TAM_ACUERDO
    print(io, COLORES, "\\begin{tikzpicture}\n\\begin{axis}[name=principal, ", EJE, ", width=", s.w, "pt, height=",
          s.h, "pt, ylabel={caudal \$U\$}, legend pos=south west]\n")
    curva(io, "k-", d.pmp.t, d.pmp.U; ancho = 1.8)
    print(io, "\\addlegendentry{disparo: PMP, 2 incógnitas}\n")
    let (t, q) = escalera(d.directo.t, d.directo.q)
        curva(io, "C0--", t, q; ancho = 1.0, extra = ", const plot")
    end
    print(io, "\\addlegendentry{directo: Ipopt, \$N=", d.directo.N, "\$}\n")
    curva(io, "C2:", d.hjb.t, d.hjb.U; ancho = 1.0, extra = ", const plot")
    print(io, "\\addlegendentry{HJB: malla \$", d.hjb.nl, "^2\\times", d.hjb.nt, "\$, lazo cerrado}\n")
    print(io, "\\end{axis}\n")
    print(io, "\\begin{axis}[at={(principal.north east)}, anchor=north east, xshift=-5pt, yshift=-5pt, ", EJE_CHICO,
          ", width=", s.wr, "pt, height=", s.hr, "pt, ylabel={\$|\\Delta U|\$}, ",
          "ylabel style={yshift=-2pt}, xtick={0,1,2}, enlarge x limits=0.03]\n")
    # Piso de 1e-16 para que un cero exacto no rompa la escala logarítmica. Ningún valor medido lo
    # alcanza; está por si una corrida futura da una coincidencia exacta.
    curva(io, "C0--", d.dif_directo.t, max.(d.dif_directo.d, 1e-16); ancho = 0.8)
    curva(io, "C2:", d.dif_hjb.t, max.(d.dif_hjb.d, 1e-16); ancho = 0.8)
    print(io, "\\end{axis}\n\\end{tikzpicture}%\n")
end

function escribir_figura_acuerdo(d, dir)
    mkpath(dir)
    archivo = joinpath(dir, "acuerdo.tex")
    open(archivo, "w") do io
        println(io, "% Escrito por el paquete Fractura (src/figuras.jl). No editar a mano.")
        figura_acuerdo(io, d)
    end
    println("figuras: acuerdo.tex, $(round(Int, filesize(archivo) / 1024)) KB")
    archivo
end

# ---- F3.3 · la frontera control ↔ información ----------------------------------------------------
# Dos paneles. El izquierdo NO grafica log det M, que es el criterio que se optimiza, sino su
# consecuencia legible: la cota de Cramér–Rao RELATIVA de cada parámetro a lo largo de la frontera,
# contra el costo de control J, que es el precio que se paga. Así se lee de una vez lo que el
# resultado dice —tres de los cuatro parámetros ya están identificados al 1–2 % en el óptimo de
# control puro, y el único que hay que comprar es η— en vez de un número (log det M) que nadie lee.
# Escala logarítmica en y porque las cotas abarcan dos órdenes y en lineal las tres chicas se
# apoyarían en el eje. Identidad por color Y por trazo, y η en negro grueso: es la protagonista.
# El panel derecho son los controles de los tres γ que se citan en el texto, y muestra el hallazgo
# que no está en los números: lo que cambia con γ es la FORMA del control (rampa decreciente →
# pulso tardío saturado), no su tamaño.
# El punto γ = 0,3 va marcado aparte, en rojo y hueco: agotó las 500 iteraciones con ‖g_proy‖ =
# 1,2e-1, así que es un iterado y no un óptimo, y la figura no lo puede presentar como uno.
# Un solo eje y por panel (nada de dos escalas en un marco), ymax con aire arriba para que la nota
# de lectura no se monte sobre la curva de η, y la leyenda en la franja vacía del medio.
const TAM_FRONTERA = (w = 152.0, h = 122.0, sep = 46.0)
# En el orden de PARS_SENS; `pars` del json lo confirma en tiempo de escritura.
const TRAZOS_PARS = ("C0--", "C2:", "k-", "C1-.")
const ANCHOS_PARS = (1.0, 1.0, 1.8, 1.0)
const ROTULOS_PARS = Dict("n" => "n", "c" => "c", "η" => "\\eta", "P0" => "P_0")
const GAMMAS_ROTULADAS = (codo = GAMMA_CODO, iterado = 0.3, extremo = 3.0)

# El punto de la frontera con ese γ. Comparación exacta: los γ vienen de GAMMAS_FRONTERA, que son
# los mismos literales, y un punto que no esté es un error de datos, no algo que redondear.
punto_gamma(fr, γ) = fr[findfirst(r -> Float64(r["gamma"]) == γ, fr)]

function curva_marcada(io, trazo, x, y; ancho = 1.2, marca = "*", tam = 1.5, extra = "")
    print(io, "\\addplot[", TRAZO[trazo], ", line width=", ancho, "pt, mark=", marca, ", mark size=", tam,
          "pt, mark options={solid}", extra, "] table[row sep=\\\\] {x y\\\\\n")
    for (a, b) in zip(x, y)
        @printf(io, "%.6g %.6g\\\\\n", a, b)
    end
    print(io, "};\n")
end

marca_suelta(io, x, y, opciones) =
    print(io, "\\addplot[", opciones, ", only marks] table[row sep=\\\\] {x y\\\\\n",
          @sprintf("%.6g %.6g\\\\\n", x, y), "};\n")

function panel_cotas(io, o)
    fr = o["frontera"]
    pars = String.(o["diseno"]["pars"])
    J = [Float64(r["J"]) for r in fr]
    cota(k) = [100 * Float64(r["cota_relativa"][k]) for r in fr]
    k_eta = findfirst(isequal("η"), pars)
    r0, rc, ri, rx = fr[1], punto_gamma(fr, GAMMAS_ROTULADAS.codo), punto_gamma(fr, GAMMAS_ROTULADAS.iterado),
                     punto_gamma(fr, GAMMAS_ROTULADAS.extremo)
    eta(r) = 100 * Float64(r["cota_relativa"][k_eta])
    print(io, "\\nextgroupplot[xlabel={costo de control \$J\$}, ylabel={cota relativa [\\%]}, ymode=log, ",
          "enlarge x limits=0.08, enlarge y limits=false, ymin=0.15, ymax=250, ytick={0.3,1,3,10,30,100}, ",
          "yticklabels={0{,}3,1,3,10,30,100}, legend style={at={(0.98,0.74)}, anchor=north east, font=\\tiny, draw=none, ",
          # ⚠️ 15-nov-2026: `fill=white, fill opacity=0.85` se PROBÓ y se revirtió. Tapa el problema
          # —la muestra negra de η se leía como un tramo de curva— pero el velo blanco DESTIÑE la
          # curva negra real y el marcador del codo que pasan por detrás: quedan grises y se leen
          # como otra serie. Peor que el defecto. El arreglo correcto es MOVER la leyenda al hueco
          # de arriba a la derecha, no ponerle fondo. Detalle en docs/marco/figuras-a-la-vista/.
          #
          # 15-nov-2026, 23.º relevo: movida. En (0.03,0.5) las CUATRO muestras quedaban sobre el
          # área de datos y las cuatro se leían como tramos de curva —no solo la negra de η—,
          # porque `curva_marcada` les pone marcador y la muestra de `c` caía encima del trazo de
          # la curva negra. Medido a 400 dpi sobre la plana 10. Se mueve al hueco de arriba a la
          # derecha (anchor north east en (0.98,0.74)), bajo el nodo de anotación y sobre las
          # curvas, que ahí no pasan de cota 5. Se mantiene `draw=none, fill=none` como todas las
          # demás leyendas del informe: si estuviera vacío el hueco, no hace falta fondo.
          "fill=none, column sep=5pt}, legend columns=2]\n")
    for k in eachindex(pars)
        curva_marcada(io, TRAZOS_PARS[k], J, cota(k); ancho = ANCHOS_PARS[k])
        print(io, "\\addlegendentry{\$", ROTULOS_PARS[pars[k]], "\$}\n")
    end
    # El codo y el iterado, sobre la curva de η. Van después de las curvas para quedar encima, y
    # fuera de la leyenda (`forget plot`).
    marca_suelta(io, Float64(rc["J"]), eta(rc), "black, mark=square, mark size=3pt, line width=0.9pt, forget plot")
    marca_suelta(io, Float64(ri["J"]), eta(ri), "rojo, mark=o, mark size=3pt, line width=0.9pt, forget plot")
    print(io, "\\node[anchor=north west, align=left, font=\\tiny, inner sep=1pt] at (rel axis cs:0.015,0.985) {",
          "\$\\gamma=0\$, a la izquierda: el óptimo de control puro\\\\\n",
          @sprintf("\$\\square\$ codo \$\\gamma=0{,}03\$: \$+%d\$\\,\\%% de \$J\$ y \$\\eta\$ %s\$\\times\$ mejor\\\\\n",
                   round(Int, 100 * (Float64(rc["J"]) / Float64(r0["J"]) - 1)), dec(eta(r0) / eta(rc), 1)),
          @sprintf("extremo \$\\gamma=3\$: \$+%d\$\\,\\%% de \$J\$ y %s\$\\times\$\\\\\n",
                   round(Int, 100 * (Float64(rx["J"]) / Float64(r0["J"]) - 1)), dec(eta(r0) / eta(rx), 1)),
          "{\\color{rojo}\$\\circ\$} \$\\gamma=0{,}3\$: iterado, no óptimo};\n")
end

function panel_controles(io, o)
    p = Par()
    fr = o["frontera"]
    print(io, "\\nextgroupplot[ylabel={caudal \$U\$}, enlarge y limits=false, ymin=-0.08, ymax=", p.qmax + 0.3,
          ", enlarge x limits=0.03, xtick={0,1,2}, legend style={at={(0.02,0.98)}, anchor=north west, ",
          "font=\\scriptsize, draw=none, fill=none}]\n")
    for (γ, trazo, rotulo, ancho) in ((0.0, "k-", "\$\\gamma=0\$: control puro", 1.6),
                                      (GAMMAS_ROTULADAS.codo, "C0--", "\$\\gamma=0{,}03\$: el codo", 1.1),
                                      (GAMMAS_ROTULADAS.extremo, "C1-.", "\$\\gamma=3\$: diseño", 1.1))
        q = Float64.(punto_gamma(fr, γ)["q"])
        let (t, u) = escalera(malla(p, length(q)), q)
            curva(io, trazo, t, u; ancho, extra = ", const plot")
        end
        print(io, "\\addlegendentry{", rotulo, "}\n")
    end
    print(io, "\\draw[gris40, dotted, line width=0.8pt] ({rel axis cs:0,0} |- {axis cs:0,", p.qmax,
          "}) -- ({rel axis cs:1,0} |- {axis cs:0,", p.qmax, "});\n")
    print(io, "\\node[anchor=north east, xshift=-2pt, yshift=-2pt, font=\\scriptsize, text=gris30] at (axis cs:",
          p.T, ",", p.qmax, ") {\$q_{\\max}\$};\n")
end

function figura_frontera(io, o)
    s = TAM_FRONTERA
    print(io, COLORES, "\\begin{tikzpicture}\n\\begin{groupplot}[group style={group size=2 by 1, horizontal sep=",
          s.sep, "pt}, ", EJE, ", width=", s.w, "pt, height=", s.h, "pt]\n")
    panel_cotas(io, o)
    panel_controles(io, o)
    print(io, "\\end{groupplot}\n\\end{tikzpicture}%\n")
end

function escribir_figura_frontera(o, dir)
    mkpath(dir)
    archivo = joinpath(dir, "frontera.tex")
    open(archivo, "w") do io
        println(io, "% Escrito por el paquete Fractura (src/figuras.jl). No editar a mano.")
        figura_frontera(io, o)
    end
    println("figuras: frontera.tex, $(round(Int, filesize(archivo) / 1024)) KB")
    archivo
end

# ---- F1.4 · el LQ de seguimiento, que es lo que justifica haber linealizado ----------------------
# Cuatro paneles: una columna por perturbación (P(0) un 20 % arriba, y η un 20 % arriba, que es el
# mismo parámetro que F3.3 muestra que NO se puede medir — por eso este experimento importa), y dos
# filas, el largo ℓ(τ) arriba y el caudal U(τ) abajo. Tres trazos por panel: la nominal (lo que se
# quería), el lazo abierto sobre el sistema perturbado (lo que pasa si uno no corrige) y el LQ de
# seguimiento (lo que pasa si corrige). La lectura está en el ESPACIO entre la nominal y el lazo
# abierto al final del horizonte, así que el panel de ℓ lleva marcado ese error terminal.
# Las series son las de `datos_figura_lq`, que rehace el experimento con `casos_lq`, la misma lista
# que usa el json: la figura y las cifras no pueden estar mirando dos experimentos distintos.
const TAM_LQ = (w = 150.0, h = 92.0, sep = 44.0, vsep = 26.0)
const ROTULOS_CASOS = Dict("P0" => "\$P(0)\$ un 20\\,\\% arriba", "eta" => "\$\\eta\$ un 20\\,\\% arriba")

function figura_lq(io, d)
    s = TAM_LQ
    print(io, COLORES, "\\begin{tikzpicture}\n\\begin{groupplot}[group style={group size=2 by 2, horizontal sep=",
          s.sep, "pt, vertical sep=", s.vsep, "pt, x descriptions at=edge bottom}, ", EJE, ", width=", s.w,
          "pt, height=", s.h, "pt]\n")
    for (fila, campo) in ((1, :ell), (2, :u))
        for (col, c) in enumerate(d.casos)
            rotulo = fila == 1 ? "largo \$\\ell\$" : "caudal \$U\$"
            extra = fila == 2 ? ", enlarge y limits=false, ymin=-0.08, ymax=$(d.qmax + 0.3)" : ""
            print(io, "\\nextgroupplot[ylabel={", rotulo, "}", extra,
                  fila == 2 && col == 1 ? ", legend pos=north west, legend style={font=\\tiny, draw=none, fill=none}" : "",
                  fila == 1 ? ", title={\\scriptsize " * ROTULOS_CASOS[c.cual] * "}" : "", "]\n")
            curva(io, "k-", d.t, getfield(d.nominal, campo); ancho = 1.6)
            fila == 2 && col == 1 && print(io, "\\addlegendentry{nominal}\n")
            curva(io, "C1-.", d.t, getfield(c.ab, campo); ancho = 1.1)
            fila == 2 && col == 1 && print(io, "\\addlegendentry{lazo abierto}\n")
            curva(io, "C0--", d.t, getfield(c.ce, campo); ancho = 1.1)
            fila == 2 && col == 1 && print(io, "\\addlegendentry{LQ de seguimiento}\n")
            if fila == 1
                # El resultado del experimento es el error TERMINAL, que en una curva de ℓ(τ) casi no
                # se ve: va escrito, con el factor que el criterio del plan exige mayor que 3.
                da = abs(c.ab.xT[1] - d.ellT_nominal)
                dc = abs(c.ce.xT[1] - d.ellT_nominal)
                print(io, "\\node[anchor=north west, align=left, font=\\tiny, inner sep=1pt] at (rel axis cs:0.03,0.97) {",
                      "desviación en \$\\tau=T\$:\\\\\n{\\color{C1}lazo abierto \$", dec(100 * da / d.ellT_nominal, 2),
                      "\$\\,\\%} · {\\color{C0}LQ \$", dec(100 * dc / d.ellT_nominal, 2), "\$\\,\\%}\\\\\n",
                      "factor \$", dec(da / dc, 1), "\$};\n")
            else
                print(io, "\\draw[gris40, dotted, line width=0.8pt] ({rel axis cs:0,0} |- {axis cs:0,", d.qmax,
                      "}) -- ({rel axis cs:1,0} |- {axis cs:0,", d.qmax, "});\n")
            end
        end
    end
    print(io, "\\end{groupplot}\n\\end{tikzpicture}%\n")
end

function escribir_figura_lq(d, dir)
    mkpath(dir)
    archivo = joinpath(dir, "lq-seguimiento.tex")
    open(archivo, "w") do io
        println(io, "% Escrito por el paquete Fractura (src/figuras.jl). No editar a mano.")
        figura_lq(io, d)
    end
    println("figuras: lq-seguimiento.tex, $(round(Int, filesize(archivo) / 1024)) KB")
    archivo
end

# ---- F4 · objetivo 5: el programa de carga óptimo y la evolución del daño ------------------------
# Dos paneles, el control y el estado: a la izquierda el programa de carga d_k que resuelve el
# seguimiento —el patrón cerrado, el descenso proyectado y el punto de partida—, a la derecha el campo
# α(x) que ese programa produce, etapa por etapa. Son el MISMO experimento: la barra con defecto en el
# centro y α = 0 en los extremos, así que el daño localiza y α(x,t) tiene forma en vez de ser una
# constante (en la barra homogénea α es uniforme y el panel derecho serían seis rectas).
# El recuadro va en el panel izquierdo porque el programa de carga es creciente y deja libre la
# esquina inferior derecha, y compara |Λ−blanco| con |Λ−blanco|: el cerrado contra el descenso, en
# escala logarítmica, que es la única forma de ver 1e-14 y 2e-4 en el mismo marco.
# No se usa groupplot: el recuadro necesita anclar a un eje con nombre.
const TAM_AT2 = (w = 148.0, h = 122.0, sep = 50.0, wr = 66.0, hr = 40.0)
# Las seis etapas, de la más clara a la más oscura: identidad por color Y por trazo.
const TRAZOS_AT2 = ("C0:", "C2:", "C0--", "C2--", "C1--", "k-")
const ANCHOS_AT2 = (0.8, 0.8, 0.9, 0.9, 1.0, 1.6)

function panel_programa_at2(io, d)
    s = TAM_AT2
    k = collect(1.0:length(d.blancos))
    print(io, "\\begin{axis}[name=programa, ", EJE, ", width=", s.w, "pt, height=", s.h,
          "pt, xlabel={etapa \$k\$}, ylabel={carga \$d_k/(\\varepsilon_c L)\$}, ",
          "xtick={1,2,3,4,5,6}, enlarge x limits=0.06, enlarge y limits=false, ymin=0, ymax=0.62, ",
          "ytick={0,0.1,0.2,0.3,0.4,0.5,0.6}, yticklabels={0,,0{,}2,,0{,}4,,0{,}6}, ",
          "legend style={at={(0.03,0.97)}, anchor=north west, font=\\scriptsize, draw=none, fill=none}]\n")
    curva(io, "C0:", k, d.carga_inicial ./ d.escala; ancho = 1.0, extra = ", mark=none")
    print(io, "\\addlegendentry{punto de partida}\n")
    curva_marcada(io, "C1--", k, d.carga_descenso ./ d.escala; ancho = 1.0, marca = "o", tam = 2.0)
    print(io, "\\addlegendentry{descenso proyectado}\n")
    curva_marcada(io, "k-", k, d.carga_cerrada ./ d.escala; ancho = 1.4, marca = "*", tam = 1.4)
    print(io, "\\addlegendentry{óptimo cerrado}\n")
    print(io, "\\end{axis}\n")
    # El recuadro: el error de seguimiento por etapa, que es lo que el gráfico grande no puede mostrar
    # porque las dos cargas se superponen.
    # `yshift=6pt` dejaba los rótulos del eje x del recuadro montados sobre la línea del eje del panel
    # grande, y `axis background` de EJE_CHICO solo pinta el ÁREA de datos: los rótulos de los ejes del
    # recuadro quedan fuera de ella y la curva punteada del punto de partida pasaba por encima del
    # «10^-8». Medido el 15-nov-2026 sobre el PDF a 400 dpi. Se sube el recuadro y se le da fondo
    # blanco a cada rótulo, que es el mismo recurso de la caja blanca de la figura del objetivo 1.
    # 17pt y no 15: con 15 la curva punteada del punto de partida asomaba 0,7pt sobre el borde
    # superior derecho del recuadro y quedaba un punto azul suelto que se leía como una mota.
    # Medido a 400 dpi el 15-nov-2026; con 17pt desaparece y el óptimo cerrado sigue a ~8pt.
    print(io, "\\begin{axis}[at={(programa.south east)}, anchor=south east, xshift=-6pt, yshift=17pt, ",
          EJE_CHICO, ", width=", s.wr, "pt, height=", s.hr,
          "pt, ylabel={\$|\\Lambda_k-\\Lambda^\\ast_k|\$}, ",
          "ylabel style={yshift=-2pt, fill=white, inner sep=1pt}, ",
          "tick label style={font=\\tiny, fill=white, inner sep=1pt}, ",
          "xtick={1,3,6}, enlarge x limits=0.08, ymin=1e-16, ymax=1e-2, ytick={1e-14,1e-8,1e-4}]\n")
    curva(io, "C1--", k, max.(abs.(d.lambda_descenso .- d.blancos), 1e-16); ancho = 0.8)
    curva(io, "k-", k, max.(abs.(d.lambda_cerrada .- d.blancos), 1e-16); ancho = 0.8)
    print(io, "\\end{axis}\n")
end

function panel_dano_at2(io, d)
    s = TAM_AT2
    print(io, "\\begin{axis}[at={(programa.south east)}, anchor=south west, xshift=", s.sep, "pt, ",
          EJE, ", width=", s.w, "pt, height=", s.h, "pt, xlabel={\$x/L\$}, ylabel={daño \$\\alpha\$}, ",
          "xtick={0,0.25,0.5,0.75,1}, xticklabels={0,,0{,}5,,1}, enlarge x limits=0.03, ",
          "ytick={0,0.05,0.1,0.15,0.2}, yticklabels={0,0{,}05,0{,}10,0{,}15,0{,}20}, ",
          "enlarge y limits=false, ymin=0, ymax=0.245, ",
          "legend style={at={(0.03,0.97)}, anchor=north west, font=\\scriptsize, draw=none, fill=none}]\n")
    x = d.x ./ d.p.L
    for k in eachindex(d.alphas)
        curva(io, TRAZOS_AT2[k], x, d.alphas[k]; ancho = ANCHOS_AT2[k])
    end
    print(io, "\\addlegendentry{\$k=1\$}\n")
    for _ in 2:(length(d.alphas) - 1)
        print(io, "\\addlegendentry{}\n")
    end
    print(io, "\\addlegendentry{\$k=", length(d.alphas), "\$}\n")
    # La flecha del tiempo, sobre el centro de la barra: el daño solo puede crecer.
    a1, aK = d.alphas[1], d.alphas[end]
    c = length(a1) ÷ 2 + 1
    print(io, "\\draw[gris40, ->, line width=0.7pt] (axis cs:0.5,", @sprintf("%.4f", a1[c] + 0.012),
          ") -- (axis cs:0.5,", @sprintf("%.4f", aK[c] - 0.012), ");\n")
    print(io, "\\node[anchor=west, font=\\scriptsize, text=gris30, inner sep=1pt] at (axis cs:0.52,",
          @sprintf("%.4f", (a1[c] + aK[c]) / 2), ") {irreversible};\n")
    # Arriba a la derecha: la campana de α cae a cero en x = L, así que esa esquina queda libre en
    # todas las etapas. Abajo a la derecha, que es donde iría por costumbre, pisa las primeras.
    print(io, "\\node[anchor=north east, align=right, font=\\scriptsize, inner sep=1pt] ",
          "at (rel axis cs:0.985,0.98) {defecto de \$G_c\$:\\\\\n",
          "\$-$(round(Int, 100 * d.p.defecto))\$\\,\\% en \$x/L=0{,}5\$,\\\\\n",
          "ancho \$", dec(d.p.ancho / d.p.L, 2), "\\,L\$};\n")
    print(io, "\\end{axis}\n")
end

function figura_at2_control(io, d)
    print(io, COLORES, "\\begin{tikzpicture}\n")
    panel_programa_at2(io, d)
    panel_dano_at2(io, d)
    print(io, "\\end{tikzpicture}%\n")
end

function escribir_figura_at2_control(d, dir)
    mkpath(dir)
    archivo = joinpath(dir, "at2-control.tex")
    open(archivo, "w") do io
        println(io, "% Escrito por el paquete Fractura (src/figuras.jl). No editar a mano.")
        figura_at2_control(io, d)
    end
    println("figuras: at2-control.tex, $(round(Int, filesize(archivo) / 1024)) KB")
    archivo
end

# ---- F0 · objetivo 1: el plano de fase, y por qué el problema de control no es trivial -----------
# Dos paneles en el MISMO plano de fase (ℓ, P), uno por familia de §3.1. Izquierda, U ≡ 0: cuatro
# presiones iniciales, cada trayectoria muriendo en su equilibrio (ℓ_∞, 0) — ahí van dos marcas, la
# raíz de la integral (3.4) y la cota explícita (3.5), de modo que la figura DIBUJA las dos
# afirmaciones derivadas en vez de ilustrarlas. Derecha, U constante: cuatro caudales hasta T,
# ninguno con equilibrio (el sistema no se detiene mientras se inyecte), y el umbral U* que alcanza
# ℓd exactamente en τ = T. La lectura de la figura es la línea ℓd, que la familia libre no cruza
# nunca y la forzada cruza solo por encima de U*: eso es lo que hace no trivial el §6.
# Las cifras van en las leyendas, no en un recuadro: son una por curva y la leyenda ya las ordena.
const TAM_OBJ1 = (w = 152.0, h = 116.0, sep = 46.0)
const TRAZOS_OBJ1 = ("C0:", "k-", "C1-.", "C2--")
const ANCHOS_OBJ1 = (1.0, 1.6, 1.1, 1.1)

"La línea vertical ℓ = ℓd: el objetivo que las dos familias cruzan o no."
linea_objetivo(io, ℓd) =
    print(io, "\\draw[gris40, dashed, line width=0.8pt] ({axis cs:", ℓd, ",0} |- {rel axis cs:0,0}) -- ",
          "({axis cs:", ℓd, ",0} |- {rel axis cs:0,1});\n")

# El rótulo ℓd va como tick extra del eje, no como nodo adentro: un nodo arriba de la línea lo tapa la
# leyenda, que pgfplots dibuja al final pase lo que pase (medido el 15-nov: la caja blanca de la leyenda
# del panel libre se comía la mitad del subíndice).
rotulo_objetivo(ℓd) = "extra x ticks={$ℓd}, extra x tick labels={\$\\ell_d\$}, " *
                      "extra x tick style={grid=none, tick style={gris40}, " *
                      "tick label style={font=\\scriptsize, text=gris30}}, "

function panel_libre_obj1(io, d)
    print(io, "\\nextgroupplot[xlabel={largo \$\\ell\$}, ylabel={presión \$P\$}, ", rotulo_objetivo(d.p.ℓd),
          "legend style={at={(0.97,0.97)}, anchor=north east, font=\\tiny, draw=none, ",
          "fill=white, fill opacity=0.85, text opacity=1}, ",
          "title={\\scriptsize sin inyectar, \$U\\equiv0\$}]\n")
    for (k, c) in enumerate(d.libres)
        curva(io, TRAZOS_OBJ1[k], c.ℓ, c.P; ancho = ANCHOS_OBJ1[k])
        print(io, "\\addlegendentry{\$P_0=", dec(c.P0, 2), "\$: \$\\ell_\\infty=", dec(c.ell_inf, 3), "\$}\n")
    end
    for c in d.libres
        marca_suelta(io, c.ell_inf, 0.0, "black, mark=*, mark size=1.5pt, forget plot")
        marca_suelta(io, c.cota, 0.0, "gris40, mark=|, mark size=3pt, line width=0.9pt, forget plot")
    end
    linea_objetivo(io, d.p.ℓd)
end

# El caudal entero se escribe entero: «U = 0,0» es ruido donde «U = 0» dice lo mismo.
rot_caudal(U) = U == round(U) ? string(round(Int, U)) : dec(U, 1)

function panel_forzada_obj1(io, d)
    print(io, "\\nextgroupplot[xlabel={largo \$\\ell\$}, ylabel={presión \$P\$}, ", rotulo_objetivo(d.p.ℓd),
          "legend style={at={(0.97,0.03)}, anchor=south east, font=\\tiny, draw=none, ",
          "fill=white, fill opacity=0.85, text opacity=1}, ",
          "title={\\scriptsize \$U\$ constante, hasta \$\\tau=T\$}]\n")
    for (k, c) in enumerate(d.forzadas)
        curva(io, TRAZOS_OBJ1[k], c.ℓ, c.P; ancho = ANCHOS_OBJ1[k])
        print(io, "\\addlegendentry{\$U=", rot_caudal(c.U), "\$: \$\\ell(T)=", dec(c.ellT, 3), "\$}\n")
        marca_suelta(io, c.ellT, c.PT, "black, mark=*, mark size=1.5pt, forget plot")
    end
    # El umbral no es un recuadro de texto: es una trayectoria más, la única que termina SOBRE la línea
    # ℓd. Dibujarla es lo que convierte «hay umbral» en algo que se ve.
    curva(io, "rojo-", d.umbral.ℓ, d.umbral.P; ancho = 1.3)
    print(io, "\\addlegendentry{\$U^*=", dec(d.Ustar, 3), "\$: \$\\ell(T)=\\ell_d\$}\n")
    marca_suelta(io, d.umbral.ellT, d.umbral.PT, "rojo, mark=*, mark size=1.5pt, forget plot")
    linea_objetivo(io, d.p.ℓd)
end

function figura_obj1(io, d)
    s = TAM_OBJ1
    print(io, COLORES, "\\begin{tikzpicture}\n\\begin{groupplot}[group style={group size=2 by 1, horizontal sep=",
          s.sep, "pt}, ", EJE, ", width=", s.w, "pt, height=", s.h, "pt]\n")
    panel_libre_obj1(io, d)
    panel_forzada_obj1(io, d)
    print(io, "\\end{groupplot}\n\\end{tikzpicture}%\n")
end

function escribir_figura_obj1(d, dir)
    mkpath(dir)
    archivo = joinpath(dir, "obj1-plano-fase.tex")
    open(archivo, "w") do io
        println(io, "% Escrito por el paquete Fractura (src/figuras.jl). No editar a mano.")
        figura_obj1(io, d)
    end
    println("figuras: obj1-plano-fase.tex, $(round(Int, filesize(archivo) / 1024)) KB")
    archivo
end
