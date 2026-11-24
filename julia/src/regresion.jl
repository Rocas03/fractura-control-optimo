# Antes de sobrescribir avance/resultados/<n>.json se compara lo recién calculado con lo que hay en disco. Si el json
# no existe no hay regresión: así se acepta un cambio de modelo, borrando el json a propósito.
const RES = joinpath(RAIZ, "avance", "resultados")
# Los resultados de la entrega final van aparte, en final/resultados; los del avance no se tocan.
const RES_FINAL = joinpath(RAIZ, "final", "resultados")

# Tolerancias: una general relativa y, por ruta del json, umbrales (residuos que solo deben ser chicos),
# absolutas y relativas más holgadas. Las rutas que calzan con EXCLUIDAS no se comparan.
const RTOL, ATOL = 1e-6, 1e-12
const EXCLUIDAS = [r"/mensaje$"]
const EXTRA_PERMITIDAS = Set(["/qmax08/fijadas", "/qmax08/dJ_fijar"])
const UMBRAL = [r"residuo_balance" => 1e-9, r"norma_residuo$" => 1e-10, r"^/disparo_penalizado/H_oscilacion$" => 1e-8,
                r"^/directo/\d+/norma_grad_inf$" => 1e-8, r"^/(pref\[\d+\]|terminal|qmax08)/H_oscilacion$" => 1e-10]
const ABSOLUTA = [r"pmp_brecha_max$" => 5e-7]
const RELATIVA = [r"^/(directo/\d+|qmax08/directo_N\d+)/H_oscilacion$" => 5e-3,
                  r"^/(directo/\d+|qmax08/directo_N\d+)/q_(min|ultimo)$" => 5e-5]

# Aplana el json en hojas con rutas «/clave» para diccionarios y «[i]», con i desde 0, para listas.
function hojas!(out, x, ruta = "")
    if x isa AbstractDict
        for (k, v) in x
            hojas!(out, v, "$ruta/$k")
        end
    elseif x isa AbstractVector
        for (i, v) in enumerate(x)
            hojas!(out, v, "$ruta[$(i - 1)]")
        end
    else
        out[ruta] = x
    end
    out
end

esnum(x) = x isa Number && !(x isa Bool)
buscar(reglas, ruta) = (i = findfirst(r -> occursin(r.first, ruta), reglas); i === nothing ? nothing : reglas[i].second)

# Lista de hojas que faltan, sobran o se salen de su tolerancia.
function fuera_de_tolerancia(ref, nuevo)
    r, n = hojas!(Dict{String, Any}(), ref), hojas!(Dict{String, Any}(), nuevo)
    malas = String[]
    for ruta in sort!(collect(union(keys(r), keys(n))))
        any(p -> occursin(p, ruta), EXCLUIDAS) && continue
        if !haskey(n, ruta)
            push!(malas, "$ruta  falta en la nueva")
            continue
        end
        if !haskey(r, ruta)
            ruta in EXTRA_PERMITIDAS || push!(malas, "$ruta  sobra en la nueva")
            continue
        end
        a, b = r[ruta], n[ruta]
        if !(esnum(a) && esnum(b))
            isequal(a, b) || push!(malas, "$ruta  distinta: $(repr(a)) vs $(repr(b))")
            continue
        end
        if (u = buscar(UMBRAL, ruta)) !== nothing
            ok, regla = abs(b) < u, "umbral $u"
        elseif (t = buscar(ABSOLUTA, ruta)) !== nothing
            ok, regla = abs(a - b) <= t, "abs $t"
        elseif (t = buscar(RELATIVA, ruta)) !== nothing
            ok, regla = abs(a - b) <= max(ATOL, t * abs(a)), "rel $t"
        else
            ok, regla = abs(a - b) <= max(ATOL, RTOL * abs(a)), "general"
        end
        ok || push!(malas, @sprintf("%s  %s: ref=%.10e nueva=%.10e", ruta, regla, a, b))
    end
    length(r), malas
end

# Se compara lo que se va a escribir (JSON.json de o, leído de vuelta), no el Dict de Julia.
function regresion(nombre, o)
    archivo = joinpath(RES, nombre * ".json")
    isfile(archivo) || return println("$nombre: no hay json en disco, no hay regresión")
    total, malas = fuera_de_tolerancia(JSON.parsefile(archivo), JSON.parse(JSON.json(o)))
    foreach(m -> println("  FUERA  ", m), malas)
    println("$nombre: $total hojas en disco, $(length(malas)) fuera de tolerancia")
    isempty(malas) || error("regresión: $nombre cambió fuera de tolerancia")
end
