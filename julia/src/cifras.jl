# Cifras del reporte de avance: un \newcommand por cada número, leído de los json recién escritos. Antes de
# escribir comprueba que lo que el texto afirma con palabras siga siendo cierto (la cota no se activa en el
# caso base, la brecha directo-disparo baja como h², todos los criterios en verde, etc.); si algo falla, da
# error y no escribe nada. Es una traducción fiel de un script anterior en Python, y por eso imita su
# formato (`pot10`, `lista_py`).
const CIFRAS = joinpath(RAIZ, "avance", "reporte", "cifras.tex")

coma(s) = replace(s, "." => "{,}")
# x con k decimales y coma decimal.
dec(x, k) = coma(Printf.format(Printf.Format("%.$(k)f"), x))
# Parámetro tal como se escribe: 2, 0,2, 10.
g(x) = coma(@sprintf("%g", x))
# Como en Python: 10**e es entero exacto si e ≥ 0 y pow() de libm si e < 0; 10.0^-2 de Julia no da lo mismo.
pot10(e) = e >= 0 ? Float64(big(10)^e) : ccall(:pow, Float64, (Float64, Float64), 10.0, e)

# x en notación a·10^e con `cifras` cifras significativas; techo = true redondea hacia arriba.
function sci(x, cifras = 1; techo = false)
    e = floor(Int, log10(abs(x)))
    m = x / pot10(e)
    esc = 10^(cifras - 1)
    m = (techo ? ceil(m * esc - 1e-9) : round(m * esc)) / esc
    if m >= 10
        m, e = m / 10, e + 1
    end
    coma(Printf.format(Printf.Format("%.$(cifras - 1)f"), m)) * "\\cdot10^{$e}"
end

pct(x, k) = dec(100 * x, k)
conjunto(valores) = "\\{" * join((g(v) for v in valores), ";\\,") * "\\}"
exigir(cond, mensaje) = cond || error("cifras_reporte: $mensaje")
lista_py(xs) = "[" * join(("'$x'" for x in xs), ", ") * "]"

function cifras_reporte(res = RES)
    leer(nombre) = JSON.parsefile(joinpath(res, nombre * ".json"))
    b1, b3, b4 = leer("objetivo1_dinamica"), leer("objetivo3_pontryagin"), leer("objetivo4_hjb")
    for (nombre, b) in (("objetivo 1", b1), ("objetivo 3", b3), ("objetivo 4", b4))
        malos = [k for (k, v) in b["criterios"] if !v]
        exigir(isempty(malos), "$nombre: criterios en rojo $(lista_py(malos))")
    end

    c = Pair{String, String}[]   # en orden: un Dict lo desordena
    poner(k, v) = push!(c, k => v)
    p = b3["parametros_base"]
    exigir(p["n"] == 2, "el texto dice «con n = 2 el campo es C^infinito»: n cambió")
    exigir(p["Pref"] == 0, "el costo (2) del reporte no tiene P_ref")
    for k in ("n", "c", "eta", "ell0", "P0", "T", "Ld", "wT", "wp", "wq", "qmax")
        poner("par" * replace(k, "0" => "cero"), g(p[k]))
    end
    d200 = b3["directo"]["200"]
    Umax = b3["disparo_penalizado"]["U_max"]
    exigir(d200["n_cota_sup"] == 0 && Umax < p["qmax"], "la cota q_max se activa en el caso base")
    exigir(d200["n_cota_inf"] == 0, "el óptimo del caso base toca U = 0 antes de T")
    poner("maxU", dec(Umax, 2))

    # Objetivo 1
    bar = b1["barrido"]
    exigir(minimum(f["min_P"] for f in bar) > 0, "P no es positiva en todo el barrido")
    exigir(minimum(f["min_diff_ell"] for f in bar) >= 0, "ell no es creciente en todo el barrido")
    poner("nSim", string(length(bar)))
    poner("nConj", conjunto(sort(unique(f["n"] for f in bar))))
    poner("etaConj", conjunto(sort(unique(f["eta"] for f in bar))))
    poner("UConj", conjunto(sort(unique(f["U"] for f in bar))))
    poner("Tbarrido", g(b1["T_barrido"]))
    poner("resBalance", sci(maximum(f["residuo_balance"] for f in bar); techo = true))

    # Objetivo 3: disparo y directo
    dp = b3["disparo_penalizado"]
    Js = dp["J"]
    poner("Jstar", dec(Js, 6))
    poner("resDisparo", sci(dp["norma_residuo"]; techo = true))
    poner("Hosc", sci(dp["H_oscilacion"]; techo = true))
    Ns = sort(collect(keys(b3["directo"])); by = N -> parse(Int, N))
    brechas = [b3["directo"][N]["J_menos_disparo"] for N in Ns]
    razones = Set(dec(brechas[i] / brechas[i + 1], 1) for i in 1:length(brechas) - 1)
    exigir(razones == Set([dec(4.0, 1)]), "la brecha directo-disparo no baja 4 veces por duplicación: $razones")
    poner("razonBrecha", only(razones))
    poner("Nprimero", Ns[1])
    poner("Nultimo", Ns[end])
    poner("brechaPrimero", sci(brechas[1], 2))
    poner("brechaUltimo", sci(brechas[end], 2))

    # Tabla
    cmp_ = b3["comparaciones"]
    filas = [
        ("Directo, \$N=$(Ns[end])\$ (aprox.~del extremal)", cmp_["optimo_directo_N200_desglose"], 6, nothing),
        ("Mejor rampa lineal escalonada, \$N=$(N_RAMPA)\$",
         cmp_["rampa"]["desglose"], 6, cmp_["rampa"]["J"]),
        ("Mejor constante, \$U=$(dec(cmp_["constante"]["U"], 3))\$",
         cmp_["constante"]["desglose"], 6, cmp_["constante"]["J"]),
        ("Constante de igual volumen, \$U=$(dec(cmp_["igual_volumen"]["U"], 3))\$",
         cmp_["igual_volumen"]["desglose"], 6, cmp_["igual_volumen"]["J"]),
    ]
    exigir(Js < cmp_["rampa"]["J"] <= cmp_["constante"]["J"] < cmp_["igual_volumen"]["J"],
           "el orden óptimo < rampa < constante < igual volumen no se cumple")
    lineas = String[]
    for (nombre, d, kJ, J) in filas
        if J === nothing
            sobre = "---"
        else
            r = (J - Js) / Js
            sobre = pct(r, r < 0.01 ? 2 : r < 0.1 ? 1 : 0) * "\\,\\%"
        end
        push!(lineas, join([nombre, dec(d["J"], kJ), sobre, dec(d["terminal"], 4),
                            dec(d["presion"], 4), dec(d["control"], 4), dec(d["volumen"], 3),
                            dec(d["ell_T"], 3)], " & ") * "\\\\")
    end
    poner("filasTabla", join(lineas, "\n"))
    poner("Ufactible", dec(cmp_["factible"]["U"], 3))

    iv = cmp_["igual_volumen"]
    poner("igualVolPct", pct((iv["J"] - Js) / Js, 0))
    poner("picoP", dec(dp["pico_P"]["P"], 3))
    poner("picoT", dec(dp["pico_P"]["t"], 2))
    exigir(cmp_["constante"]["desglose"]["P_T"] > dp["P_T"], "la mejor constante no deja más presión en T")
    poner("PTconst", dec(cmp_["constante"]["desglose"]["P_T"], 2))
    poner("PTopt", dec(dp["P_T"], 2))
    exigir(iv["desglose"]["ell_T"] < cmp_["optimo_directo_N200_desglose"]["ell_T"],
           "la constante de igual volumen no deja la grieta más corta")
    poner("ellIgual", dec(iv["desglose"]["ell_T"], 3))
    poner("rampaA", dec(cmp_["rampa"]["a"], 2))
    poner("rampaB", dec(cmp_["rampa"]["b"], 2))
    poner("rampaPct", pct((cmp_["rampa"]["J"] - Js) / Js, 2))

    # Cota activa
    q8 = b3["qmax08"]
    exigir(length(q8["conmutaciones"]) == 1, "con q_max = 0,8 no hay exactamente una conmutación")
    exigir(abs(q8["t_salida_cota_directo"] - q8["conmutaciones"][1]) <= q8["h"],
           "el directo no sale de la cota a una celda")
    poner("qmaxCota", g(q8["directo_N200"]["q_max"]))
    poner("conmutacion", dec(q8["conmutaciones"][1], 3))
    poner("dJcota", sci(abs(q8["dJ"])))
    poner("cotaPct", pct((q8["J"] - Js) / Js, 2))

    # Variante terminal
    te = b3["terminal"]
    poner("Jterminal", dec(te["J"], 6))
    poner("nu", dec(te["nu_multiplicador"], 3))
    wts = [f["wT"] for f in te["continuacion"]]
    poner("wTprimero", g(wts[1]))
    poner("wTultimo", "10^{$(round(Int, log10(wts[end])))}")
    exigir(all(8.5 <= r <= 11.5 for r in te["razones_por_decada"]), "el error ya no baja ≈10 veces por década")
    poner("razonesDecada", join((dec(r, 1) for r in te["razones_por_decada"]), ";\\ "))

    # Objetivo 4
    niv = b4["niveles"]
    dl, dP = niv[1]["dominio_ell"], niv[1]["dominio_P"]
    exigir(all(f["dominio_ell"] == dl && f["dominio_P"] == dP for f in niv), "el dominio cambia entre mallas")
    exigir(dP[2] > 1.5 * dp["pico_P"]["P"], "el dominio en P no cubre el pico de presión con margen")
    poner("domEll", "[$(g(dl[1]));\\,$(g(dl[2]))]")
    poner("domP", "[$(g(dP[1]));\\,$(g(dP[2]))]")
    ordenes = [log2(q) for q in b4["cocientes_e_V"]]
    exigir(all(0.8 <= o <= 1.3 for o in ordenes), "el orden de V_h no es cercano a 1: $ordenes")
    for (i, f) in enumerate(niv)
        poner("errV$("ABC"[i])", pct(f["e_V_rel"], 2))
        poner("malla$("ABC"[i])", "$(f["nl"])^2\\times$(f["nt"])")
        exigir(f["nl"] == f["nP"], "malla no cuadrada en (ell, P)")
    end
    poner("dJlazo", sci(niv[end]["dJ_lc"], 2))
    z = niv[end]["dzeta_max_sobre_z0"]
    length(z) == 2 || error("cifras_reporte: dzeta_max_sobre_z0 no tiene 2 componentes")
    poner("zetaEll", pct(z[1], 1))
    poner("zetaP", pct(z[2], 1))
    c
end

# Escribe una macro \cf<nombre> por cifra, comprobando que cada nombre sea solo letras.
function escribir_cifras(c, archivo = CIFRAS; etiqueta = "cifras_reporte")
    lineas = ["% Escrito por julia/src (paquete Fractura) a partir de los json de resultados. No editar a mano."]
    for (k, v) in c
        exigir(!isempty(k) && all(isletter, k), "nombre de macro inválido: $k")
        push!(lineas, "\\newcommand{\\cf$k}{$v}")
    end
    mkpath(dirname(archivo))
    write(archivo, join(lineas, "\n") * "\n")
    println("$etiqueta: $(length(c)) cifras en $(basename(archivo))")
end
