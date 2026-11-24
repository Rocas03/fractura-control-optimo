# Objetivo 5 (C5): control óptimo del programa de carga en el modelo de daño AT2, con sus criterios y las
# cifras que van al póster y al informe. Es el análogo en EDP de Pontryagin en el modelo reducido: un
# gradiente verificado contra diferencias finitas, un programa óptimo que coincide con uno cerrado y el
# punto donde la aplicación carga → estado no es diferenciable, medido como un quiebre.
# `criterios` son las condiciones que tienen que cumplirse; `lectura` es lo que se lee del resultado, que
# puede salir para cualquier lado sin que falle nada.
# Hay dos barras: la de control (defecto al centro, α = 0 en los extremos), donde el daño se localiza y
# α(x, t) tiene forma, y la homogénea, la única con solución cerrada, que sirve de control del forward.

# Los experimentos, fijados acá y no dentro de las funciones.

const P_AT2_CONTROL = ParAT2(; L = 1.0, ℓ = 0.25, N = 60, defecto = 0.4, ancho = 0.1, fija_borde = true)
const BLANCOS_AT2 = collect(range(0.02, 0.12; length = 6))          # ∫α dx pedido, etapa por etapa
const ITER_AT2 = 300
const MALLA_AT2 = (; L = 1.0, ℓ = 0.25, defecto = 0.4, ancho = 0.1, fija_borde = true)
const NS_MALLA_AT2 = (40, 80, 160, 320)
const KS_ETAPAS_AT2 = (1, 2, 3, 5, 10)              # etapas iguales hasta 0,9 ε_c L
const RAMPA_AT2 = (desde = 0.1, hasta = 1.3, etapas = 30)           # carga en unidades de ε_c L
const P_QUIEBRE_AT2 = ParAT2(; L = 1.0, ℓ = 0.3, N = 40, defecto = 0.0)
const CARGA_QUIEBRE_AT2 = [0.4, 0.6, 0.8, 0.8, 0.9]                 # la etapa 4 repite la carga de la 3
const OBJ_QUIEBRE_AT2 = (objetivo = 0.25, omega = 0.2)
const TAUS_QUIEBRE_AT2 = (1e-5, 1e-6, 1e-7)
const P_GRAD_AT2 = ParAT2(; L = 1.0, ℓ = 0.3, N = 40, defecto = 0.0)
const P_HOMOGENEA_AT2 = ParAT2(; L = 1.0, ℓ = 0.3, N = 32, defecto = 0.0)
# Rango donde la rama homogénea es la solución del esquema.
const RAMA_HOMOGENEA_AT2 = [0.2, 0.5, 0.8, 1.0, 1.4]
# Escaneo que busca el salto a daño total. El punto del salto depende del paso de carga, así que va al
# json junto con su rampa.
const ETAPAS_ESCANEO_AT2 = 47

rampa_at2(q) = collect(range(RAMPA_AT2.desde, RAMPA_AT2.hasta; length = RAMPA_AT2.etapas)) .* eps_c(q) .* q.L

# Carga de ruptura estructural: el máximo de σ en una rampa monótona hasta 1,3 ε_c L. Es lo que converge
# en malla y lo que escala como ℓ^(-1/2).
function pico_at2(q::ParAT2)
    s = resolver_at2(rampa_at2(q), q)
    k = argmax(s.sigma)
    (sigma = s.sigma[k], etapa = k, carga = s.carga[k], kkt = s.kkt, alpha_final = s.alpha[end])
end

# Convergencia en malla del pico: razones entre diferencias sucesivas y el orden que implican.
function orden_malla_at2(; ns = NS_MALLA_AT2)
    ps = [pico_at2(ParAT2(; MALLA_AT2..., N = n)) for n in ns]
    σ = [r.sigma for r in ps]
    difs = abs.(diff(σ))
    razones = difs[1:(end - 1)] ./ difs[2:end]
    (ns = collect(ns), picos = σ, difs = difs, razones = razones, ordenes = log2.(razones),
     orden = minimum(log2.(razones)), kkt = maximum(r.kkt for r in ps))
end

# σ_pico·√ℓ con L/ℓ, h/ℓ y ancho/ℓ fijos: la ley ℓ^(-1/2) del AT2, sin ajustar ninguna recta.
# `sin_raiz` es el control de mutación: sin el √ℓ los dos picos difieren claramente.
function invariancia_escala_at2()
    a = ParAT2(; L = 1.0, ℓ = 0.25, defecto = 0.4, ancho = 0.1, N = 80, fija_borde = true)
    b = ParAT2(; L = 2.0, ℓ = 0.50, defecto = 0.4, ancho = 0.2, N = 80, fija_borde = true)
    pa, pb = pico_at2(a).sigma, pico_at2(b).sigma
    (pico_a = pa, pico_b = pb, ell_a = a.ℓ, ell_b = b.ℓ,
     relativo = abs(pa * sqrt(a.ℓ) - pb * sqrt(b.ℓ)) / (pa * sqrt(a.ℓ)),
     sin_raiz = abs(pa - pb) / pa)
end

# Control analítico del forward: en la barra homogénea con cotas naturales la rama uniforme es punto
# estacionario exacto, así que el solver tiene que reproducir α(ε), σ(ε) y el pico σ_c a tolerancia de
# solver, no de discretización.
function rama_homogenea_at2(q::ParAT2 = ParAT2(; L = 1.0, ℓ = 1.0, N = 40, defecto = 0.0))
    εs = RAMA_HOMOGENEA_AT2 .* eps_c(q)
    s = resolver_at2(εs .* q.L, q)
    eα = maximum(abs(maximum(s.alpha[k]) - alpha_homogenea(εs[k], q)) for k in eachindex(εs))
    eσ = maximum(abs(s.sigma[k] - sigma_homogenea(εs[k], q)) for k in eachindex(εs))
    pk = resolver_at2([eps_pico_homogenea(q) * q.L], q)
    # La rama homogénea es punto estacionario para todo ε, pero deja de ser el mínimo que sigue el esquema:
    # con una rampa larga la barra salta a daño total y σ ya no es la de la rama. Por eso el criterio se
    # mide en el rango fijado arriba y el salto se reporta aparte.
    largo = collect(range(0.2, 2.5; length = ETAPAS_ESCANEO_AT2)) .* eps_c(q)
    sl = resolver_at2(largo .* q.L, q)
    k = findfirst(a -> maximum(a) >= q.amax - 1e-9, sl.alpha)
    ult = findlast(j -> abs(maximum(sl.alpha[j]) - alpha_homogenea(largo[j], q)) < 1e-9,
                   eachindex(largo))
    (error_alpha = eα, error_sigma = eσ,
     error_pico = abs(pk.sigma[end] - sigma_pico_homogenea(q)) / sigma_pico_homogenea(q),
     kkt = max(s.kkt, pk.kkt), rango = collect(RAMA_HOMOGENEA_AT2),
     salto = k === nothing ? NaN : largo[k] / eps_c(q),
     ultima_que_calza = ult === nothing ? NaN : largo[ult] / eps_c(q))
end

# Gradiente adjunto contra diferencias finitas centradas en tres direcciones, donde el conjunto activo no
# cambia. `mutacion` es el gradiente sin el término de rango uno del Hessiano reducido, que sale con el
# signo cambiado.
function gradiente_contra_fd_at2()
    p = P_GRAD_AT2
    εc = eps_c(p)
    K = 8
    blancos = collect(range(0.06, 0.40; length = K))
    om = 0.01
    carga = collect(range(0.2, 0.95; length = K)) .* εc .* p.L
    r = gradiente_at2(carga, p; objetivo = blancos, omega = om)
    costo(c) = costo_at2(resolver_at2(c, p); objetivo = blancos, omega = om).J
    τ = 1e-6 * εc * p.L
    dirs = (normalize(ones(K)), normalize(collect(1.0:K)), [k == 4 ? 1.0 : 0.0 for k in 1:K])
    errores = Float64[]
    for dir in dirs
        fd = (costo(carga .+ τ .* dir) - costo(carga .- τ .* dir)) / (2τ)
        push!(errores, abs(dot(r.grad, dir) - fd) / max(1.0, abs(fd)))
    end
    rs = gradiente_at2(carga, p; objetivo = blancos, omega = om, rango1 = false)
    dir = normalize(ones(K))
    fd = (costo(carga .+ τ .* dir) - costo(carga .- τ .* dir)) / (2τ)
    (errores = errores, congelados = r.congelados, J = r.J,
     mutacion = (adjunto = dot(rs.grad, dir), fd = fd, mismo_signo = dot(rs.grad, dir) * fd > 0,
                 error = abs(dot(rs.grad, dir) - fd) / max(1.0, abs(fd))))
end

# El quiebre: en una etapa que repite la carga anterior, α queda en la cota de irreversibilidad con
# multiplicador cero, y las derivadas laterales del costo respecto de esa carga difieren en un factor que
# no se mueve al bajar τ cien veces. No es error de truncamiento: ahí la aplicación carga → estado no es
# diferenciable, el adjunto da la derivada de la rama y el principio del máximo clásico no está disponible.
function quiebre_at2()
    p = P_QUIEBRE_AT2
    εc = eps_c(p)
    carga = CARGA_QUIEBRE_AT2 .* εc .* p.L
    obj, om = OBJ_QUIEBRE_AT2.objetivo, OBJ_QUIEBRE_AT2.omega
    r = gradiente_at2(carga, p; objetivo = obj, omega = om)
    costo(c) = costo_at2(resolver_at2(c, p); objetivo = obj, omega = om).J
    e4 = [k == 4 ? 1.0 : 0.0 for k in eachindex(carga)]
    base = costo(carga)
    taus = Dict{String, Any}[]
    for f in TAUS_QUIEBRE_AT2
        τ = f * εc * p.L
        der = (costo(carga .+ τ .* e4) - base) / τ
        izq = (base - costo(carga .- τ .* e4)) / τ
        push!(taus, Dict("tau_relativo" => f, "derecha" => der, "izquierda" => izq,
                         "razon" => izq / der, "centrada" => (der + izq) / 2))
    end
    razones = [t["razon"] for t in taus]
    # Control: con la etapa 4 un poco más arriba el conjunto activo queda vacío y el mismo código calza con
    # diferencias finitas. El desacuerdo es del punto, no del adjunto.
    limpia = [0.4, 0.6, 0.8, 0.85, 0.9] .* εc .* p.L
    r2 = gradiente_at2(limpia, p; objetivo = obj, omega = om)
    τ = 1e-6 * εc * p.L
    fd2 = (costo(limpia .+ τ .* e4) - costo(limpia .- τ .* e4)) / (2τ)
    # Separación del adjunto a cada lateral: a diferencia de la distancia al promedio, no depende de cómo el
    # redondeo desempata el conjunto activo.
    sep = min(abs(r.grad[4] / taus[end]["derecha"]), abs(taus[end]["izquierda"] / r.grad[4]))
    (taus = taus, razones = razones, adjunto = r.grad[4], congelados = r.congelados,
     deriva_de_la_razon = maximum(razones) / minimum(razones) - 1,
     promedio = taus[end]["centrada"], separacion_laterales = sep,
     control_limpio = (congelados = r2.congelados, adjunto = r2.grad[4], fd = fd2,
                       error = abs(r2.grad[4] - fd2) / max(1.0, abs(fd2))))
end

# El problema de control: llevar ∫α dx al programa pedido etapa por etapa. Con omega = 0 el óptimo es
# cerrado (invertir Λ(d), que es creciente), así que el descenso proyectado se mide contra una referencia
# conocida y no contra sí mismo.
function seguimiento_at2(p::ParAT2, blancos::AbstractVector; iter = ITER_AT2)
    K = length(blancos)
    dex = carga_exacta_at2(blancos, p)
    sx = resolver_at2(dex, p)
    v0 = fill(dex[end] / K * 0.6, K)
    o = optimizar_at2(v0, p; objetivo = blancos, omega = 0.0, iter)
    g = gradiente_at2(dex, p; objetivo = blancos, omega = 0.0)
    (blancos = collect(blancos), carga_cerrada = dex, carga_descenso = o.carga,
     carga_inicial = cumsum(v0), lambda_cerrada = sx.longitud, lambda_descenso = o.sol.longitud,
     sigma_cerrada = sx.sigma, alphas = sx.alpha, x = nodos_at2(p), escala = eps_c(p) * p.L,
     error_cerrado = maximum(abs, sx.longitud .- blancos),
     # `error_cerrado` es la precisión de la bisección que arma la referencia; esta es la distancia del
     # descenso a ella, programa a programa.
     distancia_al_cerrado_pct = 100 * maximum(abs.(o.carga .- dex) ./ dex),
     J_cerrado = costo_at2(sx; objetivo = blancos, omega = 0.0).J,
     error_descenso = maximum(abs, o.sol.longitud .- blancos),
     relativo_descenso = maximum(abs, o.sol.longitud .- blancos) / blancos[end],
     pasos = length(o.historia), J0 = o.historia[1], J = o.J, historia = o.historia,
     costo_baja = all(diff(o.historia) .< 0), carga_monotona = all(diff(o.carga) .>= 0),
     congelados = g.congelados, fijos = count(nodos_fijos(p)), kkt = sx.kkt, p = p)
end

# Descenso con blanco terminal (escalar) hacia el mismo blanco final del seguimiento, con omega = 0. Ahí el
# descenso sí alcanza el blanco: lo que el blanco terminal deja indeterminado es el programa, no el punto
# final, y por eso el costo es de seguimiento.
function deficit_terminal_at2(p::ParAT2, blancos::AbstractVector)
    K = length(blancos)
    term = blancos[end]
    dex = carga_exacta_at2(blancos, p)
    o = optimizar_at2(fill(dex[end] / K * 0.6, K), p; objetivo = term, omega = 0.0, iter = ITER_AT2)
    lam = o.sol.longitud[end]
    (blanco = term, longitud = lam, deficit_pct = 100 * (1 - lam / term),
     pasos = length(o.historia))
end

# Residuo KKT de llegar a la misma carga final (0,9 ε_c L) en la barra `P_AT2_CONTROL`, repartida en K
# etapas iguales. Muestra que la minimización alternada necesita pasos de carga chicos: con dos y con tres
# etapas el residuo todavía queda sobre 1e-2.
function kkt_por_etapas_at2(; ks = KS_ETAPAS_AT2)
    p = P_AT2_CONTROL
    dfin = 0.9 * eps_c(p) * p.L
    Dict(string(K) => resolver_at2(collect(range(dfin / K, dfin; length = K)), p).kkt for K in ks)
end

# El objetivo 5 completo: criterios, lectura y las series de la figura. Devuelve `(o, ct)`: `o` va al json
# y `ct` a la figura de `escribir_figura_at2_control`, igual que `acuerdo()`.
function objetivo5_at2()
    hom = rama_homogenea_at2()
    # La rama homogénea por defecto tiene ℓ = L, y el punto del salto depende de ℓ/L; esta segunda barra es
    # la de la figura y la del experimento de malla.
    hom_fig = rama_homogenea_at2(ParAT2(; L = 1.0, ℓ = MALLA_AT2.ℓ, N = 40, defecto = 0.0))
    ml = orden_malla_at2()
    es = invariancia_escala_at2()
    gr = gradiente_contra_fd_at2()
    qb = quiebre_at2()
    ct = seguimiento_at2(P_AT2_CONTROL, BLANCOS_AT2)
    hg = seguimiento_at2(P_HOMOGENEA_AT2, collect(range(0.06, 0.38; length = 6)))
    tc = deficit_terminal_at2(P_AT2_CONTROL, BLANCOS_AT2)
    th = deficit_terminal_at2(P_HOMOGENEA_AT2, collect(range(0.06, 0.38; length = 6)))
    kkt = kkt_por_etapas_at2()

    criterios = Dict(
        # Forward: el control analítico y la convergencia en malla.
        "la_rama_homogenea_reproduce_la_cerrada_bajo_1e9" =>
            max(hom.error_alpha, hom.error_sigma, hom.error_pico) < 1e-9,
        "el_orden_de_convergencia_en_malla_es_al_menos_085" => ml.orden >= 0.85,
        "la_invariancia_de_escala_pasa_1e9" => es.relativo < 1e-9,
        "mutacion_sin_el_raiz_de_ell_los_picos_difieren_mas_de_20pct" => es.sin_raiz > 0.2,
        # Adjunto: el gradiente contra diferencias finitas y su mutación.
        "el_gradiente_adjunto_concuerda_con_diferencias_finitas_bajo_1e5" =>
            maximum(gr.errores) < 1e-5,
        "mutacion_sin_el_termino_de_rango_uno_el_gradiente_cambia_de_signo" =>
            !gr.mutacion.mismo_signo,
        # El quiebre, que es el resultado.
        "las_derivadas_laterales_difieren_mas_de_3x" => minimum(qb.razones) > 3,
        "la_razon_de_las_laterales_no_se_mueve_1pct_al_bajar_tau_100x" =>
            qb.deriva_de_la_razon < 0.01,
        # El adjunto tiene que quedar separado de las dos laterales por al menos un 25 %; compararlo con el
        # promedio de las laterales era frágil frente al redondeo.
        "el_adjunto_no_es_ninguna_de_las_dos_laterales" => qb.separacion_laterales > 1.25,
        "donde_el_conjunto_activo_no_cambia_el_mismo_codigo_concuerda" =>
            qb.control_limpio.error < 1e-5 && all(==(0), qb.control_limpio.congelados),
        # Control óptimo: el patrón cerrado y el descenso, en las dos barras.
        "el_optimo_cerrado_alcanza_el_blanco_bajo_1e10" =>
            max(ct.error_cerrado, hg.error_cerrado) < 1e-10,
        "el_descenso_queda_bajo_el_5pct_del_blanco_en_las_dos_barras" =>
            max(ct.relativo_descenso, hg.relativo_descenso) < 0.05,
        "el_costo_baja_monotono_y_la_carga_sale_monotona" =>
            ct.costo_baja && ct.carga_monotona && hg.costo_baja && hg.carga_monotona,
        "el_conjunto_activo_del_optimo_solo_tiene_los_nodos_pinchados" =>
            all(==(ct.fijos), ct.congelados) && all(==(0), hg.congelados),
        "el_residuo_kkt_de_todos_los_forward_queda_bajo_1e10" =>
            max(hom.kkt, ml.kkt, ct.kkt, hg.kkt) < 1e-10)

    lectura = Dict(
        "orden_observado_en_malla" => ml.orden,
        "razones_de_las_diferencias" => ml.razones,
        "factor_entre_las_derivadas_laterales" => qb.razones[1],
        "adjunto_en_el_quiebre" => qb.adjunto,
        "promedio_de_las_laterales" => qb.promedio,
        "adjunto_relativo_al_promedio" => qb.adjunto / qb.promedio,
        # Inestable frente a perturbaciones del tamaño del redondeo: se registra, no se cita.
        "separacion_del_adjunto_a_la_lateral_mas_cercana" => qb.separacion_laterales,
        "nodos_congelados_en_la_etapa_que_sostiene" => qb.congelados[4],
        "nodos_de_la_barra_del_quiebre" => P_QUIEBRE_AT2.N + 1,
        # Déficit respecto del blanco final cuando el blanco es terminal.
        "deficit_pct_con_blanco_terminal_con_defecto" => tc.deficit_pct,
        "deficit_pct_con_blanco_terminal_homogenea" => th.deficit_pct,
        "pasos_del_descenso_con_blanco_terminal" => [tc.pasos, th.pasos],
        "error_del_seguimiento_pct_con_defecto" => 100 * ct.relativo_descenso,
        "error_del_seguimiento_pct_homogenea" => 100 * hg.relativo_descenso,
        "factor_de_bajada_del_costo_con_defecto" => ct.J0 / ct.J,
        "gradiente_error_relativo_maximo" => maximum(gr.errores),
        "pico_de_ruptura_malla_fina" => ml.picos[end],
        "pico_homogeneo_mismo_ell" => sigma_pico_homogenea(ParAT2(; ℓ = MALLA_AT2.ℓ)),
        "sin_el_raiz_de_ell_los_picos_difieren_pct" => 100 * es.sin_raiz,
        # Límite de la minimización alternada, no criterio: misma barra y misma carga final repartida en K
        # etapas. Lo que manda es el tamaño del paso de carga.
        "kkt_si_una_sola_etapa_salta_a_09_epsc" => kkt[string(1)],
        "kkt_misma_barra_dos_etapas" => kkt[string(2)],
        "kkt_misma_barra_tres_etapas" => kkt[string(3)],
        "kkt_misma_barra_cinco_etapas" => kkt[string(5)],
        "kkt_misma_barra_diez_etapas" => kkt[string(10)],
        # Ojo: este es el máximo sobre las cuatro mallas de la rampa de 30 etapas hasta 1,3 ε_c, con otra
        # barra y otra carga final que las de arriba.
        "kkt_de_la_rampa_de_30_etapas_peor_malla" => ml.kkt,
        "distancia_del_descenso_al_programa_cerrado_pct_con_defecto" => ct.distancia_al_cerrado_pct,
        "distancia_del_descenso_al_programa_cerrado_pct_homogenea" => hg.distancia_al_cerrado_pct,
        # Hasta dónde vale la rama homogénea como solución, no como punto estacionario.
        "la_rama_homogenea_calza_hasta_epsc" => hom.ultima_que_calza,
        "salta_a_dano_total_en_epsc" => hom.salto)

    o = Dict(
        "malla" => Dict("N" => ml.ns, "picos" => ml.picos, "diferencias" => ml.difs,
                        "razones" => ml.razones, "ordenes" => ml.ordenes,
                        "ell" => MALLA_AT2.ℓ, "defecto" => MALLA_AT2.defecto,
                        "ancho" => MALLA_AT2.ancho),
        "escala" => Dict("pico_a" => es.pico_a, "pico_b" => es.pico_b, "ell_a" => es.ell_a,
                         "ell_b" => es.ell_b, "relativo" => es.relativo, "sin_raiz" => es.sin_raiz),
        "homogenea" => Dict("error_alpha" => hom.error_alpha, "error_sigma" => hom.error_sigma,
                            "error_pico" => hom.error_pico, "constante_del_pico" => AT2_PICO,
                            "rango_en_epsc" => hom.rango,
                            "salta_a_dano_total_en_epsc" => hom.salto,
                            "ultima_eps_que_calza_en_epsc" => hom.ultima_que_calza,
                            "barra" => Dict("L" => 1.0, "ell" => 1.0, "N" => 40),
                            "con_ell_de_la_figura" =>
                                Dict("ell" => MALLA_AT2.ℓ, "N" => 40,
                                     "salta_a_dano_total_en_epsc" => hom_fig.salto,
                                     "ultima_eps_que_calza_en_epsc" => hom_fig.ultima_que_calza),
                            "escaneo" => Dict("desde" => 0.2, "hasta" => 2.5,
                                              "etapas" => ETAPAS_ESCANEO_AT2,
                                              "unidad" => "epsilon_c")),
        "gradiente" => Dict("errores_relativos" => gr.errores, "congelados" => gr.congelados,
                            "J" => gr.J,
                            "mutacion_sin_rango_uno" =>
                                Dict("adjunto" => gr.mutacion.adjunto, "fd" => gr.mutacion.fd,
                                     "mismo_signo" => gr.mutacion.mismo_signo)),
        "quiebre" => Dict("taus" => qb.taus, "adjunto" => qb.adjunto,
                          "congelados" => qb.congelados,
                          "deriva_de_la_razon" => qb.deriva_de_la_razon,
                          "control_limpio" =>
                              Dict("congelados" => qb.control_limpio.congelados,
                                   "adjunto" => qb.control_limpio.adjunto,
                                   "fd" => qb.control_limpio.fd,
                                   "error_relativo" => qb.control_limpio.error),
                          "carga" => CARGA_QUIEBRE_AT2, "objetivo" => OBJ_QUIEBRE_AT2.objetivo,
                          "omega" => OBJ_QUIEBRE_AT2.omega),
        "control" => Dict(nombre =>
                              Dict("blancos" => s.blancos, "carga_cerrada" => s.carga_cerrada,
                                   "carga_descenso" => s.carga_descenso,
                                   "lambda_cerrada" => s.lambda_cerrada,
                                   "lambda_descenso" => s.lambda_descenso,
                                   "error_cerrado" => s.error_cerrado, "J_cerrado" => s.J_cerrado,
                                   "error_descenso" => s.error_descenso,
                                   "relativo_descenso" => s.relativo_descenso,
                                   "J0" => s.J0, "J" => s.J, "pasos" => s.pasos,
                                   "congelados" => s.congelados, "kkt" => s.kkt,
                                   "distancia_al_cerrado_pct" => s.distancia_al_cerrado_pct,
                                   "barra" => Dict("L" => s.p.L, "ell" => s.p.ℓ, "N" => s.p.N,
                                                   "defecto" => s.p.defecto, "ancho" => s.p.ancho,
                                                   "fija_borde" => s.p.fija_borde))
                          for (nombre, s) in ("con_defecto" => ct, "homogenea" => hg)),
        "kkt_por_etapas" => kkt,
        "criterios" => criterios, "lectura" => lectura,
        "modelo" => Dict("energia" => "AT2: w(α) = α², degradación (1-α)², c_w = 2",
                         "forward" => "minimización alternada; paso en u cerrado (σ constante), " *
                                      "paso en α QP convexo con cotas por conjunto activo",
                         "irreversibilidad" => "cota inferior nodo por nodo, α ≥ α de la etapa anterior",
                         "adjunto" => "discreto, hacia atrás, conjunto activo congelado; Hessiano " *
                                      "reducido con el término de rango uno de la eliminación de u",
                         "rampa" => Dict("desde" => RAMPA_AT2.desde, "hasta" => RAMPA_AT2.hasta,
                                         "etapas" => RAMPA_AT2.etapas,
                                         "unidad" => "epsilon_c * L")))
    (o, ct)
end
