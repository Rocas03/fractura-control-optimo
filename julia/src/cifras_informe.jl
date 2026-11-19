# Las cifras del INFORME FINAL: un \newcommand por cada número de los cuatro json de `final/resultados/`.
# Mismo patrón que `cifras.jl` (que hace las del avance) y reusa su formato: `coma`, `dec`, `sci`, `exigir`,
# `escribir_cifras`. La diferencia es el criterio de qué entra: aquí NO se eligen cifras a mano, se recorren
# los json completos. Así ninguna cifra del informe se tipea a mano, y si un json cambia, el macro cambia solo.
#
# Qué entra y qué no:
#   · todo valor ESCALAR NUMÉRICO (Int o Float), recorriendo los diccionarios anidados;
#   · los vectores (y las listas de diccionarios) de a lo más MAX_COMPONENTES componentes, componente por
#     componente —esto incluye las matrices 2×2, que son vectores de vectores—;
#   · NO entran: la sección "criterios" (booleanos), cualquier otro booleano, los strings, y los vectores de
#     más de MAX_COMPONENTES componentes (series: `nominal/q`, `frontera`, `barrido_condicionamiento`, …).
#     `cifras_informe(; detalle = true)` devuelve, junto a las cifras, la lista de todo lo que quedó fuera
#     con su razón, que es lo que se reporta al revisar la cobertura.
#
# REGLA CLAVE→NOMBRE (mecánica, documentada acá porque el nombre del macro es la única interfaz con el .tex):
# nombre = \cf + prefijo del json + un trozo CamelCase por cada clave de la ruta; cada clave se corta por todo
# lo que no sea [A-Za-z0-9], cada trozo va con mayúscula inicial (sin tocar las demás letras: `N400` sigue con
# su N) y cada DÍGITO se escribe con su palabra en mayúscula inicial (0→Cero … 9→Nueve), porque \newcommand no
# acepta dígitos ni guiones. La i-ésima componente de un vector es el trozo `cI`, o sea `CUno`, `CDos`, …
# Ejemplos: acuerdo_tres_metodos/acuerdo_en_J/J_directo_N400 → \cfacuAcuerdoEnJJDirectoNCuatroCeroCero;
#           objetivo2_lineal/gramiano_controlabilidad/lambda[1] → \cfdosGramianoControlabilidadLambdaCUno.
# La unicidad del nombre se EXIGE: dos rutas distintas que den el mismo nombre son error, no silencio.
const CIFRAS_FINAL = joinpath(RAIZ, "final", "informe", "cifras-final.tex")

# Los cinco json de la entrega final, cada uno con su prefijo de macro. El quinto, `cruce_lyapunov_atribucion`,
# lo escribe julia/cruce_atribucion.jl (19-nov): la atribución del cruce de Lyapunov del §4.4, que antes vivía
# tipeada a mano desde docs/marco/ataque-s5/.
const JSON_INFORME = (("objetivo2_lineal", "dos"), ("acuerdo_tres_metodos", "acu"),
                      ("diseno_fisher", "fis"), ("objetivo5_at2", "atf"), ("cruce_lyapunov_atribucion", "cru"))
# El informe final cita además UN bloque que NO está en final/resultados/: la prueba del borde artificial de
# HJB (hallazgo H14 de la revisión de Astra), medida en avance/resultados/objetivo4_hjb.json y nunca publicada.
# `avance/` es entrega cerrada: acá solo se LEE. Prefijo propio `hjb` para no chocar con cifras-avance.tex.
const RES_AVANCE = joinpath(RAIZ, "avance", "resultados")
const JSON_AVANCE_INFORME = (("objetivo4_hjb", "hjb", ["borde"]),)
const MAX_COMPONENTES = 4
const PALABRA_DIGITO = ("Cero", "Uno", "Dos", "Tres", "Cuatro", "Cinco", "Seis", "Siete", "Ocho", "Nueve")

# Un trozo de la ruta, en CamelCase ASCII y sin dígitos.
function trozo_nombre(clave::AbstractString)
    partes = split(clave, r"[^A-Za-z0-9]+"; keepempty = false)
    sin_digitos(p) = join(isdigit(ch) ? PALABRA_DIGITO[Int(ch) - Int('0') + 1] : ch for ch in p)
    join(uppercasefirst(sin_digitos(p)) for p in partes)
end
nombre_cifra(prefijo, ruta) = prefijo * join(trozo_nombre(k) for k in ruta)

# El número tal como va al informe: entero tal cual; float con 4 cifras significativas, en notación científica
# si es muy grande o muy chico (los residuos del orden de 1e-15 no se escriben con coma decimal).
function numero_informe(x)
    x isa Integer && return string(x)
    exigir(isfinite(x), "valor no finito: $x")
    x == 0 && return "0"
    a = abs(x)
    (a >= 1e4 || a < 1e-3) && return sci(x, 4)
    isinteger(x) && return g(x)
    dec(x, max(0, 3 - floor(Int, log10(a))))
end

function _recorrer!(c, fuera, prefijo, ruta, v)
    nombre() = nombre_cifra(prefijo, ruta)
    camino() = join(ruta, "/")
    if v isa AbstractDict
        for k in sort(collect(keys(v)))
            k == "criterios" && (push!(fuera, (camino() * "/" * k, "sección de criterios (booleanos)")); continue)
            _recorrer!(c, fuera, prefijo, [ruta; k], v[k])
        end
    elseif v isa AbstractVector
        if length(v) <= MAX_COMPONENTES
            for (i, x) in enumerate(v)
                _recorrer!(c, fuera, prefijo, [ruta; "c$i"], x)
            end
        else
            push!(fuera, (camino(), "vector de $(length(v)) componentes (> $MAX_COMPONENTES)"))
        end
    elseif v isa Bool
        push!(fuera, (camino(), "booleano"))
    elseif v isa Number
        push!(c, nombre() => numero_informe(v))
    elseif v isa AbstractString
        push!(fuera, (camino(), "string"))
    else
        push!(fuera, (camino(), "tipo no numérico: $(typeof(v))"))
    end
end

function cifras_informe(res = RES_FINAL; detalle = false)
    c = Pair{String, String}[]   # en orden: un Dict lo desordena
    fuera = Tuple{String, String}[]
    for (archivo, prefijo) in JSON_INFORME
        b = JSON.parsefile(joinpath(res, archivo * ".json"))
        exigir(b isa AbstractDict, "$archivo.json no es un diccionario")
        _recorrer!(c, fuera, prefijo, String[], b)
    end
    for (archivo, prefijo, ruta) in JSON_AVANCE_INFORME
        b = JSON.parsefile(joinpath(RES_AVANCE, archivo * ".json"))
        exigir(b isa AbstractDict, "$archivo.json no es un diccionario")
        sub = b
        for k in ruta
            exigir(sub isa AbstractDict && haskey(sub, k),
                   "cifras_informe: $archivo.json no tiene la ruta $(join(ruta, "/"))")
            sub = sub[k]
        end
        _recorrer!(c, fuera, prefijo, String[ruta...], sub)
    end
    visto = Dict{String, Int}()
    for (i, (k, _)) in enumerate(c)
        exigir(occursin(r"^[A-Za-z]+$", k), "nombre de macro inválido (solo letras ASCII): $k")
        # El mensaje se arma aparte: `exigir` evalúa sus dos argumentos, y visto[k] todavía no existe.
        haskey(visto, k) && error("cifras_informe: nombre de macro repetido: \\cf$k (ya salió en la cifra $(visto[k]))")
        visto[k] = i
    end
    detalle ? (c, fuera) : c
end
