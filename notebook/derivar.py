#!/usr/bin/env python3
"""Arma notebook/fractura.ipynb, el notebook autocontenido que se entrega, a partir del paquete Fractura (julia/).

El código vive una sola vez, en el paquete: cada celda de código es, tal cual, un archivo de julia/src/ o de
julia/test/, o un bloque «# %% <id>» de julia/regenerar.jl, y las versiones de la celda Entorno salen del [compat] de
julia/Project.toml. ESQUELETO fija el orden de las celdas y los títulos. Sobre ese notebook base van los cambios de la
entrega.

Uso, desde cualquier carpeta:  python3 -I notebook/derivar.py [salida.ipynb]
Por defecto escribe notebook/fractura.ipynb sin salidas. Se ejecuta aparte, en una carpeta vacía nueva (jupyter
nbconvert --to notebook --execute --inplace), y ese notebook ejecutado, con sus salidas y figuras, reemplaza al
de notebook/: es el que se entrega.

Diferencias con el paquete, y nada más:
- en vez del entorno del paquete, la celda Entorno fija las mismas versiones en un entorno temporal;
- todo lo que escribe va a salidas/, junto al notebook (nunca a la raíz del repositorio, RAIZ);
- las tres figuras del reporte también se dibujan con Plots + GR, en pantalla y en salidas/figuras/*.png;
- la última celda compara las cifras recién calculadas con las del reporte (avance/reporte/cifras.tex, copiadas
  aquí) y da error si alguna difiere.

Cada cambio es un reemplazo anclado por contenido: la celda se ubica por un texto que aparece en una sola celda, y el
texto reemplazado tiene que aparecer exactamente una vez en ella. Si no, aborta sin escribir nada.
"""
import hashlib
import tomllib
import json
import re
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
JL = RAIZ / "julia"
CIFRAS = RAIZ / "avance" / "reporte" / "cifras.tex"
SALIDA = Path(sys.argv[1]) if len(sys.argv) > 1 else RAIZ / "notebook" / "fractura.ipynb"


PAQUETES = ("OrdinaryDiffEqVerner", "NLPModels", "NLPModelsIpopt", "QuadGK", "JSON", "StaticArrays")

ENTORNO = """\
# Versiones exactas en un entorno temporal: sin red si la caché de paquetes alcanza, con red si no.
import Pkg
const VERSIONES = [{versiones}]
Pkg.activate(; temp = true, io = devnull)
let specs = [Pkg.PackageSpec(; name, version) for (name, version) in VERSIONES]
    try
        Pkg.offline(true)
        Pkg.add(specs; io = devnull)
    catch
        Pkg.offline(false)
        Pkg.add(specs; io = devnull)
    end
end
using OrdinaryDiffEqVerner, NLPModels, NLPModelsIpopt, QuadGK, JSON, StaticArrays,
      LinearAlgebra, Printf, Random, Test"""

# (id, contenido): un .jl es una celda de código con ese archivo; «corrida», el bloque de julia/regenerar.jl con ese
# id; «entorno», ENTORNO; cualquier otra cosa es el texto de una celda markdown.
ESQUELETO = [
    ("c001", "# Fractura hidráulica: control óptimo del caudal de inyección"),
    ("c002", "## Entorno"), ("c003", "entorno"),
    ("c004", "## Parámetros y modelo"), ("c005", "src/modelo.jl"),
    ("c006", "## Integrador y lazo abierto"), ("c007", "src/integrador.jl"),
    ("c038", "## Raíces"), ("c039", "src/raices.jl"),
    ("c008", "## Método directo: avance, adjunto e Ipopt"), ("c009", "src/directo.jl"),
    ("c010", "## Disparo sobre el sistema extremal"), ("c011", "src/disparo.jl"),
    ("c012", "## Comparaciones con controles simples"), ("c013", "src/comparaciones.jl"),
    ("c014", "## HJB por semi-Lagrangiano"), ("c015", "src/hjb.jl"),
    ("c016", "## Objetivo 1: dinámica en lazo abierto"), ("c017", "src/objetivo1.jl"),
    ("c018", "## Objetivo 3: Pontryagin contra el directo y contra controles simples"), ("c019", "src/objetivo3.jl"),
    ("c020", "## Objetivo 4: HJB contra el disparo"), ("c021", "src/objetivo4.jl"),
    ("c050", "## Objetivo 2: el sistema linealizado y el LQ de seguimiento"), ("c051", "src/objetivo2.jl"),
    ("c052", "## El acuerdo entre los tres métodos"), ("c053", "src/acuerdo.jl"),
    ("c054", "## Sensibilidades paramétricas, calibre e identificabilidad"),
    ("c055", "src/sensibilidad.jl"),
    ("c056", "## Matriz de información de Fisher y cota de Cramér–Rao"),
    ("c057", "src/fisher.jl"),
    ("c058", "## La frontera control ↔ información"),
    ("c059", "src/diseno.jl"),
    ("c060", "## Objetivo 5: el daño por gradiente AT2 en 1D, el modelo de EDP"),
    ("c061", "src/at2.jl"),
    ("c062", "## El adjunto del AT2 y el gradiente respecto del programa de carga"),
    ("c063", "src/at2_adjunto.jl"),
    ("c064", "## Objetivo 5: el resultado, con sus criterios"),
    ("c065", "src/objetivo5.jl"),
    ("c022", "## Corridas"), ("c023", "corrida"), ("c024", "corrida"), ("c025", "corrida"),
    ("c040", "## Semilla: comprobaciones de lectura (numéricas en BigFloat, no simbólicas)"),
    ("c041", "src/semilla.jl"), ("c042", "corrida"),
    ("c031", "## Pruebas"), ("c032", "test/ayudas.jl"), ("c033", "test/pruebas.jl"),
    ("c034", "## Regresión contra los json en disco"), ("c035", "src/regresion.jl"), ("c036", "corrida"),
    ("c037", "corrida"),
    ("c043", "## Salidas: json, `cifras.tex` y figuras"), ("c044", "corrida"), ("c045", "src/cifras.jl"),
    ("c046", "corrida"), ("c047", "src/figuras.jl"), ("c049", "src/datos_beamer.jl"), ("c048", "corrida"),
]


def abortar(mensaje):
    sys.exit(f"derivar.py: {mensaje}; no se escribió nada")


def texto(c):
    return "".join(c["source"])


def poner(c, s):
    c["source"] = s.splitlines(keepends=True)


def ubicar(nb, ancla):
    i = [k for k, c in enumerate(nb["cells"]) if ancla in texto(c)]
    if len(i) != 1:
        abortar(f"el ancla de celda {ancla!r} aparece en {len(i)} celdas")
    return i[0]


def reemplazar(nb, ancla, viejo, nuevo):
    c = nb["cells"][ubicar(nb, ancla)]
    n = texto(c).count(viejo)
    if n != 1:
        abortar(f"{viejo[:60]!r} aparece {n} veces en la celda de {ancla!r}")
    poner(c, texto(c).replace(viejo, nuevo))


def celda(id_, tipo, s):
    c = {"cell_type": tipo, "id": id_, "metadata": {}, "source": []}
    if tipo == "code":
        c.update(execution_count=None, outputs=[])
    poner(c, s)
    return c


def entorno():
    compat = tomllib.loads((JL / "Project.toml").read_text(encoding="utf-8"))["compat"]
    if any(not compat.get(k, "").startswith("=") for k in PAQUETES):
        abortar("julia/Project.toml no fija con «=» la versión de cada paquete de PAQUETES")
    v = [f'"{k}" => "{compat[k][1:]}"' for k in PAQUETES]
    return ENTORNO.format(versiones=", ".join(v[:3]) + ",\n                   " + ", ".join(v[3:]))


def corridas():
    """Los bloques «# %% <id>» de julia/regenerar.jl, sin los saltos de línea del final."""
    partes = re.split(r"^# %% (c\d+)\n", (JL / "regenerar.jl").read_text(encoding="utf-8"), flags=re.M)
    return {partes[k]: partes[k + 1].rstrip("\n") for k in range(1, len(partes), 2)}


def base():
    """El notebook armado desde el paquete, antes de los cambios de la entrega."""
    bloques = corridas()
    usados = [i for i, x in ESQUELETO if x == "corrida"]
    if sorted(usados) != sorted(bloques):
        abortar(f"los bloques de regenerar.jl {sorted(bloques)} no son los del esqueleto {sorted(usados)}")
    nb = {"cells": [], "metadata": {"kernelspec": {"display_name": "Julia 1.11", "language": "julia",
                                                    "name": "julia-1.11"},
                                     "language_info": {"file_extension": ".jl", "mimetype": "application/julia",
                                                       "name": "julia", "version": "1.11.4"}},
          "nbformat": 4, "nbformat_minor": 5}
    for id_, x in ESQUELETO:
        if x.endswith(".jl"):
            s = (JL / x).read_text(encoding="utf-8")
            nb["cells"].append(celda(id_, "code", s[:-1] if s.endswith("\n") else s))
        elif x == "corrida":
            nb["cells"].append(celda(id_, "code", bloques[id_]))
        elif x == "entorno":
            nb["cells"].append(celda(id_, "code", entorno()))
        else:
            nb["cells"].append(celda(id_, "markdown", x))
    return nb


def insertar_despues(nb, ancla, *nuevas):
    i = ubicar(nb, ancla)
    nb["cells"][i + 1:i + 1] = list(nuevas)


def jl(s):
    """Literal de string de Julia: el escape de JSON sirve salvo por $, que en Julia interpola."""
    return json.dumps(s, ensure_ascii=False).replace("$", "\\$")


def referencia():
    """Las cifras de cifras.tex como pares (nombre, valor), con ida y vuelta exacta al texto."""
    t = CIFRAS.read_text(encoding="utf-8")
    cabecera, _, cuerpo = t.partition("\n")
    if not cabecera.startswith("% Escrito por"):
        abortar(f"cabecera inesperada en {CIFRAS}")
    pares = []
    for p in re.split(r"\n(?=\\newcommand\{\\cf)", cuerpo.rstrip("\n")):
        m = re.fullmatch(r"\\newcommand\{\\cf([A-Za-z]+)\}\{(.*)\}", p, re.S)
        if not m:
            abortar(f"no se entiende la cifra {p[:60]!r}")
        pares.append((m[1], m[2]))
    if "\n".join(f"\\newcommand{{\\cf{k}}}{{{v}}}" for k, v in pares) + "\n" != cuerpo:
        abortar("las cifras leídas no reconstruyen cifras.tex")
    return pares, hashlib.sha256(t.encode("utf-8")).hexdigest()


INTRO = """\
# Fractura hidráulica: control óptimo del caudal de inyección

Notebook autocontenido del avance del proyecto (MA4703 Control Óptimo, primavera 2026). Calcula todo lo que reporta el
avance —dinámica en lazo abierto, Pontryagin contra el método directo y HJB—, corre sus pruebas y deja los resultados en
`salidas/`, en la misma carpeta que este archivo. No necesita ningún otro archivo.

## Cómo correrlo
1. **Julia 1.11**, con juliaup: en macOS o Linux, `curl -fsSL https://install.julialang.org | sh`; en Windows, lo que
   indica https://julialang.org/install/. Después, `juliaup add 1.11` y `juliaup default 1.11`.
2. **IJulia**, el kernel de Jupyter para Julia: `julia -e 'import Pkg; Pkg.add("IJulia")'`.
3. Abrir este notebook en Jupyter (`julia -e 'using IJulia; notebook(dir = pwd())'`, desde esta carpeta) o en VS Code
   con la extensión de Julia, elegir el kernel de Julia 1.11 y **correr todas las celdas en orden**.

La celda «Entorno» fija la versión de cada paquete en un entorno temporal. La primera vez los baja de internet y los
precompila, lo que puede tardar varios minutos; las corridas siguientes no usan la red.

## Qué deja en `salidas/`
- `resultados/`: los json de resultados de los objetivos 1, 3 y 4.
- `semilla/resultados-lectura.json`: las comprobaciones de lectura del modelo.
- `reporte/cifras.tex` y `reporte/figuras/*.tex`: las cifras y las figuras que compila el reporte.
- `figuras/*.png`: las mismas tres figuras, que además se muestran en el notebook.
- `beamer/datos/*.csv`: las curvas que anima la presentación.

La última celda compara las cifras recién calculadas con las del reporte entregado y da error si alguna difiere. Tres
cifras son ruido de redondeo (residuos del orden de $10^{-12}$ a $10^{-15}$) y cambian en la última cifra según el
computador: para esas basta el orden de magnitud, y la celda lo avisa con `RUIDO` sin dar error.
"""

GRAFICOS = """\
# Las tres figuras del reporte, en pantalla (Plots + GR) y en PNG en salidas/figuras/. Mismos datos que las pgfplots
# de salidas/reporte/figuras/, que son las que compila el reporte; solo cambia quién dibuja.
const FIGURAS_PNG = joinpath(SALIDAS, "figuras")
const PALETA = Plots.palette(:tab10)   # C0, C1, C2 y el rojo de matplotlib
const GRIS = Plots.RGB(0.4, 0.4, 0.4)

function grafico_comparacion(d)
    trazos = ((:black, :solid), (PALETA[1], :dash), (PALETA[2], :dashdot), (PALETA[3], :dot))
    paneles = map(enumerate(("caudal U", "largo ℓ", "presión P"))) do (j, rotulo)
        p = Plots.plot(; xlabel = "τ", ylabel = rotulo, legend = j == 2 ? :topleft : false)
        for (c, (color, estilo)) in zip(d, trazos)
            Plots.plot!(p, c[2], c[2 + j]; label = c[1], color, linestyle = estilo, linewidth = 1.5)
        end
        p
    end
    Plots.plot(paneles...; layout = (1, 3), size = (1050, 320), left_margin = 5 * Plots.mm, bottom_margin = 5 * Plots.mm)
end

function grafico_cota(d)
    p = Plots.plot(d.t, d.U; label = "disparo", color = :black, linewidth = 1.5, xlabel = "τ", ylabel = "caudal U",
                   ylims = (-0.05, d.qmax + 0.25), legend = :bottomleft, size = (650, 340))
    Plots.plot!(p, d.td, d.Ud; label = "directo, N = $(d.N)", color = PALETA[2], linestyle = :dash, linewidth = 1.5)
    Plots.hline!(p, [d.qmax]; label = "q_max", color = GRIS, linestyle = :dot)
    Plots.vline!(p, d.conmutaciones; label = "conmutación", color = PALETA[4], linestyle = :dash)
    p
end

function grafico_retroalimentacion(d)
    p = Plots.plot(d.t, d.U; label = "U* del PMP", color = :black, linewidth = 2, xlabel = "τ", ylabel = "caudal U",
                   legend = :bottomleft, size = (650, 340))
    for (n, (color, estilo)) in zip(d.niveles, ((PALETA[1], :dot), (PALETA[2], :dashdot), (PALETA[3], :dash)))
        Plots.plot!(p, n.t, n.U; seriestype = :steppost, label = "HJB, malla $(n.nl)²×$(n.nt)", color,
                    linestyle = estilo, linewidth = 1.2)
    end
    p
end

function mostrar_figuras(d, dir = FIGURAS_PNG)
    mkpath(dir)
    for (nombre, g, x) in (("comparacion", grafico_comparacion, d.comparacion), ("cota-activa", grafico_cota, d.cota),
                           ("retroalimentacion", grafico_retroalimentacion, d.retro))
        p = g(x)
        archivo = joinpath(dir, nombre * ".png")
        Plots.savefig(p, archivo)
        display(MIME("image/png"), p)   # solo PNG: el SVG y el HTML de Plots triplican el peso
        println("figuras: $nombre.png, $(round(Int, filesize(archivo) / 1024)) KB")
    end
end
"""

COMPARAR = r"""
# Cifras que son ruido de redondeo (residuos y oscilaciones del orden de 1e-12 a 1e-15): su última cifra depende del
# procesador y de las versiones de las dependencias, no del modelo. Para ellas basta el orden de magnitud: la nueva no
# puede pasar de 10 veces la del reporte. Las demás se comparan como texto, exactas.
const RUIDO = ("resBalance", "resDisparo", "Hosc")

function valor_potencia(s)
    m = match(r"^(?:([0-9]+(?:\{,\}[0-9]+)?)\\cdot)?10\^\{(-?[0-9]+)\}$", s)
    m === nothing && error("referencia: «$s» no es de la forma a\\cdot10^{e}")
    a = m[1] === nothing ? 1.0 : parse(Float64, replace(m[1], "{,}" => "."))
    return a * 10.0^parse(Int, m[2])
end

function comparar_referencia(c, ref = REFERENCIA)
    nuevas = Dict(c)
    malas = String[]
    for (k, v) in ref
        if !haskey(nuevas, k)
            push!(malas, "$k: falta en lo recién calculado")
        elseif nuevas[k] == v
            continue
        elseif k in RUIDO && valor_potencia(nuevas[k]) <= 10 * valor_potencia(v)
            println("  RUIDO      $k: el reporte dice «$v», ahora sale «$(nuevas[k])» (mismo orden o menor: se acepta)")
        else
            push!(malas, "$k: el reporte dice «$v», ahora sale «$(nuevas[k])»")
        end
    end
    foreach(k -> push!(malas, "$k: no está en el reporte"), setdiff(first.(c), first.(ref)))
    foreach(m -> println("  DISTINTA  ", m), malas)
    println("referencia: $(length(ref)) cifras del reporte, $(length(malas)) distintas")
    isempty(malas) || error("referencia: lo calculado no coincide con el reporte entregado")
end
"""


def main():
    nb = base()
    if len(sys.argv) > 2 and sys.argv[2] == "--base":
        SALIDA.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        return print(f"derivar.py: notebook base, {len(nb['cells'])} celdas")
    pares, sha = referencia()

    poner(nb["cells"][ubicar(nb, "# Fractura hidráulica: control óptimo del caudal de inyección")], INTRO)

    reemplazar(nb, "const VERSIONES", '"StaticArrays" => "1.9.22"]',
               '"StaticArrays" => "1.9.22",\n                   "Plots" => "1.41.7", "GR" => "0.73.27"]')
    reemplazar(nb, "const VERSIONES", "      LinearAlgebra, Printf, Random, Test",
               "      LinearAlgebra, Printf, Random, Test\n"
               'ENV["GKSwstype"] = "100"   # GR sin ventana: las figuras van al notebook y a PNG\n'
               "import Plots\n\n"
               "# Todo lo que escribe el notebook va a salidas/, en la carpeta del notebook.\n"
               'const SALIDAS = joinpath(pwd(), "salidas")\n'
               "mkpath(SALIDAS);")

    reemplazar(nb, 'sem["identidades_numericas_bigfloat"]',   # único de la celda de corrida c042:
               'open(joinpath(RAIZ, "semilla", "resultados-lectura.json"), "w") do io',
               'mkpath(joinpath(SALIDAS, "semilla"))\n'
               'open(joinpath(SALIDAS, "semilla", "resultados-lectura.json"), "w") do io')

    reemplazar(nb, "function regresion(nombre, o)",
               "# Antes de sobrescribir avance/resultados/<n>.json se compara lo recién calculado con lo que hay en disco",
               "# Antes de sobrescribir salidas/resultados/<n>.json se compara lo recién calculado con lo que hay en disco")
    reemplazar(nb, "function regresion(nombre, o)",
               "Si el json no existe no hay regresión: así se\n# acepta un cambio de modelo, borrando el json a propósito.",
               "En la primera corrida no hay json y la regresión\n"
               "# queda OMITIDA; la comparación contra el reporte entregado es la de la última celda.")
    reemplazar(nb, "function regresion(nombre, o)",
               'const RES = joinpath(RAIZ, "avance", "resultados")',
               'const RES = joinpath(SALIDAS, "resultados")\nmkpath(RES)')
    reemplazar(nb, "function regresion(nombre, o)",
               '# `final/` es la entrega final: json, figuras e informe nuevos. `avance/` está entregado y no se toca.\n'
               'const RES_FINAL = joinpath(RAIZ, "final", "resultados")',
               '# En el repositorio esto es final/resultados (la entrega final); en el notebook todo va a salidas/.\n'
               'const RES_FINAL = joinpath(SALIDAS, "resultados-final")')
    reemplazar(nb, "function regresion(nombre, o)",
               'isfile(archivo) || return println("$nombre: no hay json en disco, no hay regresión")',
               'isfile(archivo) || return println("$nombre: regresión OMITIDA, no hay json de una corrida anterior")')

    reemplazar(nb, "function cifras_reporte(res = RES)",
               'const CIFRAS = joinpath(RAIZ, "avance", "reporte", "cifras.tex")',
               'const CIFRAS = joinpath(SALIDAS, "reporte", "cifras.tex")')
    reemplazar(nb, "function cifras_reporte(res = RES)", "% Escrito por julia/src (paquete Fractura)", "% Escrito por fractura.ipynb")
    reemplazar(nb, "escribir_cifras(cifras_reporte())", "escribir_cifras(cifras_reporte())",
               "cf = cifras_reporte()\nescribir_cifras(cf)")

    reemplazar(nb, "function escribir_figuras(d, dir = FIGURAS)",
               'const FIGURAS = joinpath(RAIZ, "avance", "reporte", "figuras")',
               'const FIGURAS = joinpath(SALIDAS, "reporte", "figuras")')
    reemplazar(nb, "function escribir_figuras(d, dir = FIGURAS)", "% Escrito por julia/src (paquete Fractura)",
               "% Escrito por fractura.ipynb")
    reemplazar(nb, "function escribir_datos_beamer",
               'const DATOS_BEAMER = joinpath(RAIZ, "avance", "video", "beamer", "datos")',
               'const DATOS_BEAMER = joinpath(SALIDAS, "beamer", "datos")')

    filas = "\n".join(f"    {jl(k)} => {jl(v)}," for k, v in pares)
    tabla = (f"# Las {len(pares)} cifras del reporte entregado, copiadas por notebook/derivar.py de avance/reporte/cifras.tex\n"
             f"# (sha256 {sha}). Se comparan como texto con las recién calculadas, salvo las de RUIDO (abajo),\n"
             "# que se comparan en orden de magnitud. Cualquier otra diferencia, aunque sea en la última cifra, se lista y da error. Va al final para que, aun con error, todas las\n"
             "# salidas ya estén escritas.\n"
             f"const REFERENCIA = [\n{filas}\n]\n\n" + COMPARAR)
    insertar_despues(nb, "df = datos_figuras(f3, f4, o4)",
                     celda("entrega-figuras", "markdown", "## Figuras en pantalla y en `salidas/figuras/`"),
                     celda("entrega-graficos", "code", GRAFICOS),
                     celda("entrega-mostrar", "code", "mostrar_figuras(df)"),
                     celda("entrega-referencia", "markdown", "## Referencia: las cifras del reporte entregado"),
                     celda("entrega-tabla", "code", tabla),
                     celda("entrega-comparar", "code", "comparar_referencia(cf)"))

    for c in nb["cells"]:
        c["metadata"] = {}
        if c["cell_type"] == "code":
            c["execution_count"], c["outputs"] = None, []

    listados = {x for _, x in ESQUELETO if x.startswith("src/")}
    enel = {f"src/{f.name}" for f in (JL / "src").glob("*.jl")} - {"src/Fractura.jl"}
    if enel - listados:
        abortar(f"módulos de src/ fuera del ESQUELETO: {sorted(enel - listados)}")
    if listados - enel:
        abortar(f"el ESQUELETO lista módulos que no existen: {sorted(listados - enel)}")
    sin_id = sum("id" not in c for c in nb["cells"])
    if sin_id:
        abortar(f"{sin_id} celdas sin id")
    sobra = [texto(c) for c in nb["cells"] if any(x in texto(c) for x in ('"..",', "RAIZ", "julia/", "pwd(), \"..\""))]
    if sobra:
        abortar(f"quedan rutas al repositorio en {len(sobra)} celdas")
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"derivar.py: {SALIDA.name}, {len(nb['cells'])} celdas, {len(pares)} cifras de referencia (sha256 {sha[:12]})")


main()
