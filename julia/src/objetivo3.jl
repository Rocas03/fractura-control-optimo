# Objetivo 3: el extremal de Pontryagin (disparo) contra el método directo y contra controles simples.
# SEMILLA es el ζ(0) con que parte el disparo.
const SEMILLA = [0.8161, 0.4346]
const NS = (25, 50, 100, 200)
const WTS = (10.0, 1e2, 1e3, 1e4)
const PREFS = (0.0, 0.25, 0.5)

# Método directo con N = 25, 50, 100 y 200 tramos; cada N parte del óptimo del anterior, con cada tramo partido en dos.
function cascada_directo(p)
    out, q = Dict{Int, Any}(), nothing
    for N in NS
        r = optimizar(p, N; q0 = q === nothing ? zeros(N) : repeat(q, inner = 2))
        q = r.x
        out[N] = r
    end
    out
end

sin_sol(r) = Dict{String, Any}(k => v for (k, v) in r if k != "sol")

# Curvas (t, U, ℓ, P) para las figuras: la del extremal y la de un control constante a trozos.
function curva_extremal(r, p; n = 801)
    t = range(0.0, p.T, length = n)
    Y = r["sol"].(t)
    collect(t), [ustar(x2(y), ζ2(y), p) for y in Y], getindex.(Y, 1), getindex.(Y, 2)
end

function curva_escalonada(q, p; n = 801)
    m = max(9, n ÷ length(q))
    t = malla(p, length(q))
    tramos = avance(q, p)[2]
    ts = [range(t[k], t[k+1], length = m) for k in eachindex(q)]
    Y = [tramos[k].(ts[k]) for k in eachindex(q)]
    reduce(vcat, collect.(ts)), reduce(vcat, [fill(q[k], m) for k in eachindex(q)]),
        reduce(vcat, [getindex.(y, 1) for y in Y]), reduce(vcat, [getindex.(y, 2) for y in Y])
end

# Máximo de la presión y del caudal a lo largo del extremal, en una malla fina.
function pico_presion(r, p; n = 20001)
    t = range(0.0, p.T, length = n)
    P = [r["sol"](τ)[2] for τ in t]
    i = argmax(P)
    Dict{String, Any}("t" => t[i], "P" => P[i], "n_malla" => n)
end

maximo_control(r, p; n = 20001) = maximum(curva_extremal(r, p; n)[2])

function objetivo3()
    p = Par()
    criterios = Dict{String, Any}()
    res = Dict{String, Any}("parametros_base" => comodict(p), "semilla" => SEMILLA)

    # Disparo contra directo: mismo J, el disparo bajo el directo en todo N, y ζ(0) del adjunto discreto
    # cerca del ζ(0) del disparo.
    rp = resolver(p, SEMILLA, :penalizado)
    directo = cascada_directo(p)
    q200 = directo[200].x
    z0_directo = x2(adjunto(q200, p, avance(q200, p)[2])[2][1](0.0))
    res["disparo_penalizado"] = merge(sin_sol(rp), Dict("pico_P" => pico_presion(rp, p), "U_max" => maximo_control(rp, p)))
    res["directo"] = Dict{String, Any}(string(N) => merge(Dict{String, Any}("J" => r.fun, "J_menos_disparo" => r.fun - rp["J"]),
                                                          diagnosticos(r.x, p)) for (N, r) in directo)
    res["directo"]["200"]["z0"] = z0_directo
    criterios["disparo_converge"] = rp["exito"] && rp["norma_residuo"] < 1e-10
    criterios["disparo_J_0.37186200"] = abs(rp["J"] - 0.37186200) ≤ 1e-8
    criterios["disparo_bajo_directo_todo_N"] = all(rp["J"] < r.fun for r in values(directo))
    criterios["z0_disparo_cerca_de_z0_directo_N200_5e-3"] = maximum(abs, rp["z0"] .- z0_directo) ≤ 5e-3
    criterios["H_oscilacion_1e-8"] = rp["H_oscilacion"] ≤ 1e-8

    # Controles simples: el extremal tiene que quedar bajo la mejor rampa, y esta bajo la mejor constante.
    rc, rr, rv = mejor_constante(p), mejor_rampa(p), constante_igual_volumen(p, SEMILLA)
    U_fact, ℓ_fact = constante_factible(p)
    q_rampa = base_rampa(p) * [rr["a"], rr["b"]]
    res["comparaciones"] = Dict{String, Any}(
        "optimo_disparo" => Dict("J" => rp["J"]),
        "optimo_directo_N200_desglose" => desglose(q200, p),
        "constante" => merge(rc, Dict("desglose" => desglose(fill(rc["U"], 200), p))),
        "igual_volumen" => merge(rv, Dict("desglose" => desglose(fill(rv["U"], 200), p))),
        "rampa" => merge(rr, Dict("desglose" => desglose(q_rampa, p))),
        "factible" => Dict("U" => U_fact, "ell_T" => ℓ_fact))
    criterios["factible_brentq"] = 0.0 ≤ U_fact ≤ p.qmax && abs(ℓ_fact - p.ℓd) < 1e-10
    criterios["orden_opt_rampa_const"] = rp["J"] < rr["J"] ≤ rc["J"] + 1e-12
    criterios["opt_bajo_igual_volumen"] = rp["J"] < rv["J"]
    curvas = [("extremal", curva_extremal(rp, p)), ("mejor constante", curva_escalonada([rc["U"]], p)),
              ("igual volumen", curva_escalonada([rv["U"]], p)), ("mejor rampa", curva_escalonada(q_rampa, p))]

    # Con ℓ(T) = ℓd fijo, ζℓ(T) = ν es un multiplicador libre. Al subir wT en el penalizado, -wT(ℓ(T) - ℓd)
    # debe acercarse a ν con error del orden de 1/wT, o sea unas diez veces menos por década.
    rt = resolver(p, SEMILLA, :terminal)
    z0, cont = SEMILLA, Dict{String, Any}[]
    for wT in WTS
        r = resolver(con(p; wT), z0)
        z0 = r["z0"]
        μ = -wT * (r["ell_T"] - p.ℓd)
        push!(cont, Dict{String, Any}("wT" => wT, "exito" => r["exito"], "J" => r["J"], "ell_T" => r["ell_T"],
                                      "multiplicador" => μ, "error_vs_nu" => abs(μ - rt["zeta_ell_T"])))
    end
    razones = [cont[i]["error_vs_nu"] / cont[i+1]["error_vs_nu"] for i in 1:length(cont)-1]
    directo_wT = Dict{String, Any}()
    for wT in WTS[1:2]
        pw = con(p; wT)
        r200 = cascada_directo(pw)[200]
        ℓT = diagnosticos(r200.x, pw)["ell_T"]
        directo_wT[string(wT)] = Dict{String, Any}("J" => r200.fun, "ell_T" => ℓT,
                                                   "error_vs_nu" => abs(-wT * (ℓT - p.ℓd) - rt["zeta_ell_T"]))
    end
    res["terminal"] = merge(sin_sol(rt), Dict{String, Any}("nu_multiplicador" => rt["zeta_ell_T"], "continuacion" => cont,
                                                           "razones_por_decada" => razones, "directo_N200" => directo_wT))
    criterios["terminal_ell_T_1e-10"] = rt["exito"] && abs(rt["ell_T"] - p.ℓd) < 1e-10
    criterios["continuacion_converge"] = all(c["exito"] for c in cont)
    criterios["continuacion_10x_por_decada"] = all(5 ≤ q ≤ 20 for q in razones)

    # Presión de referencia Pref distinta de cero, por continuación desde la semilla.
    z0, pref = SEMILLA, Dict{String, Any}[]
    for Pref in PREFS
        r = resolver(con(p; Pref), z0)
        z0 = r["z0"]
        push!(pref, merge(Dict{String, Any}("Pref" => Pref),
                          Dict{String, Any}(k => r[k] for k in ("exito", "J", "ell_T", "P_T", "z0", "norma_residuo", "H_oscilacion"))))
    end
    res["pref"] = pref
    criterios["pref_disparo_converge"] = all(r["exito"] for r in pref)

    # Con qmax = 0.8 la cota se activa al inicio: un solo arco en la cota, y el instante en que el directo sale
    # de ella tiene que caer a menos de un tramo h de la conmutación del disparo.
    p08 = con(p; qmax = 0.8)
    r08 = resolver(p08, SEMILLA)
    d08 = cascada_directo(p08)[200]
    q08 = d08.x
    h = p08.T / 200
    en_cota = findall(≥(p08.qmax - 1e-12), q08)
    t_directo = isempty(en_cota) ? nothing : en_cota[end] * h
    res["qmax08"] = merge(sin_sol(r08), Dict{String, Any}("directo_N200" => diagnosticos(q08, p08), "J_directo_N200" => d08.fun,
                                                          "dJ" => d08.fun - r08["J"], "t_salida_cota_directo" => t_directo,
                                                          "h" => h, "fijadas" => d08.fijadas,
                                                          "dJ_fijar" => d08.fun - d08.fun_sin_fijar))
    criterios["qmax08_converge"] = r08["exito"] && r08["norma_residuo"] < 1e-10
    criterios["qmax08_dJ_1e-5"] = abs(d08.fun - r08["J"]) ≤ 1e-5
    criterios["qmax08_un_arco_conmutacion_h"] = !isempty(en_cota) && all(diff(en_cota) .== 1) &&
                                                !isempty(r08["conmutaciones"]) && abs(r08["conmutaciones"][1] - t_directo) ≤ h
    res["criterios"] = criterios
    res, (; rp, curvas, r08, p08, q08)
end
