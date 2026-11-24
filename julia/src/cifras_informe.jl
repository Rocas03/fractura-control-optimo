# Cifras del informe final: un \newcommand por cada número de los json de `final/resultados/`, con el
# formato de `cifras.jl`. Acá no se eligen cifras a mano: se recorren los json completos, así que ninguna
# cifra del informe se tipea y, si un json cambia, la macro cambia sola.
# Entran todos los escalares numéricos y los vectores de a lo más MAX_COMPONENTES componentes, componente
# por componente (incluidas las matrices 2×2). Quedan fuera los criterios, los booleanos, los strings y las
# series largas; `cifras_informe(; detalle = true)` devuelve además lo que quedó fuera y por qué.
# El nombre de la macro es \cf + prefijo del json + un trozo CamelCase por clave de la ruta, con cada dígito
# escrito en palabra (0 → Cero, …, 9 → Nueve), porque \newcommand no acepta dígitos; la componente i de un
# vector es el trozo `cI`. Por ejemplo, acuerdo_tres_metodos/acuerdo_en_J/J_directo_N400 da
# \cfacuAcuerdoEnJJDirectoNCuatroCeroCero. Dos rutas con el mismo nombre son error.
const CIFRAS_FINAL = joinpath(RAIZ, "final", "informe", "cifras-final.tex")

# Los json de la entrega final, cada uno con su prefijo de macro; `cruce_lyapunov_atribucion` lo escribe
# julia/cruce_atribucion.jl.
const JSON_INFORME = (("objetivo2_lineal", "dos"), ("acuerdo_tres_metodos", "acu"),
                      ("diseno_fisher", "fis"), ("objetivo5_at2", "atf"), ("cruce_lyapunov_atribucion", "cru"))
# El informe cita además un bloque de avance/resultados/objetivo4_hjb.json, la prueba del borde artificial
# de HJB. `avance/` es una entrega cerrada y acá solo se lee; el prefijo `hjb` evita chocar con las cifras
# del avance.
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

# Recorre el json y separa lo que entra como cifra de lo que queda fuera, con su razón.
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
