# F4 · Objetivo 5, el resultado publicable: el control óptimo del programa de carga en el modelo de
# daño por gradiente AT2, con sus criterios y las cifras que entran al póster y al informe.
#
# Es el análogo en EDP del objetivo 3 (Pontryagin en el modelo reducido), y la razón por la que el
# objetivo 5 dejó de ser «se trasladaría»: acá hay un gradiente verificado contra diferencias finitas,
# un programa de carga óptimo que coincide con un patrón cerrado, y —el resultado que vale— el punto
# donde la aplicación carga→estado NO es diferenciable, medido como un quiebre y no como truncamiento.
#
# Los criterios son los del plan y de las pruebas 4.1/4.2, sin retocar: orden de convergencia ≥ 0,85,
# gradiente a 1e-5 relativo, seguimiento bajo el 5 % del blanco. Lo que el plan NO pedía —el factor
# entre las derivadas laterales, el orden observado, la invariancia de escala— va en «lectura», que es
# lo que se lee del resultado y puede salir para cualquier lado sin que el criterio falle.
#
# Dos barras, y la distinción importa:
#  · la de CONTROL (defecto en el centro, α = 0 en los extremos) es la del experimento de control y la
#    de la figura: el daño localiza en el defecto, así que α(x,t) tiene forma y no es una constante;
#  · la HOMOGÉNEA es la única con solución cerrada, así que es el control analítico del forward y la
#    que reproduce el 0,024 % de la prueba 4.2.

# ---- los experimentos, fijados acá y no en el cuerpo de las funciones -----------------------------

const P_AT2_CONTROL = ParAT2(; L = 1.0, ℓ = 0.25, N = 60, defecto = 0.4, ancho = 0.1, fija_borde = true)
const BLANCOS_AT2 = collect(range(0.02, 0.12; length = 6))          # ∫α dx pedido, etapa por etapa
const ITER_AT2 = 300
const MALLA_AT2 = (; L = 1.0, ℓ = 0.25, defecto = 0.4, ancho = 0.1, fija_borde = true)
const NS_MALLA_AT2 = (40, 80, 160, 320)
const KS_ETAPAS_AT2 = (1, 2, 3, 5, 10)              # etapas iguales hasta 0,9 ε_c L (B3 del ataque)
const RAMPA_AT2 = (desde = 0.1, hasta = 1.3, etapas = 30)           # carga en unidades de ε_c L
const P_QUIEBRE_AT2 = ParAT2(; L = 1.0, ℓ = 0.3, N = 40, defecto = 0.0)
const CARGA_QUIEBRE_AT2 = [0.4, 0.6, 0.8, 0.8, 0.9]                 # la etapa 4 SOSTIENE la 3
const OBJ_QUIEBRE_AT2 = (objetivo = 0.25, omega = 0.2)
const TAUS_QUIEBRE_AT2 = (1e-5, 1e-6, 1e-7)
const P_GRAD_AT2 = ParAT2(; L = 1.0, ℓ = 0.3, N = 40, defecto = 0.0)
const P_HOMOGENEA_AT2 = ParAT2(; L = 1.0, ℓ = 0.3, N = 32, defecto = 0.0)
# El rango de la prueba 4.1, sin retocar: donde la rama homogénea ES la solución del esquema.
const RAMA_HOMOGENEA_AT2 = [0.2, 0.5, 0.8, 1.0, 1.4]
# El escaneo que busca el salto. El punto del salto DEPENDE del paso de carga —con 25 etapas salta
# antes de 2,5 ε_c y con 47 aguanta hasta 2,35—, así que el número va al json junto con su rampa.
const ETAPAS_ESCANEO_AT2 = 47

rampa_at2(q) = collect(range(RAMPA_AT2.desde, RAMPA_AT2.hasta; length = RAMPA_AT2.etapas)) .* eps_c(q) .* q.L

"""
    pico_at2(q)

La carga de ruptura estructural: el máximo de σ sobre una rampa monótona hasta 1,3 ε_c L. Es la
cantidad que converge en malla y la que escala como ℓ^(-1/2).
"""
function pico_at2(q::ParAT2)
    s = resolver_at2(rampa_at2(q), q)
    k = argmax(s.sigma)
    (sigma = s.sigma[k], etapa = k, carga = s.carga[k], kkt = s.kkt, alpha_final = s.alpha[end])
end

"Convergencia en malla del pico: razones de las diferencias sucesivas y el orden que implican."
function orden_malla_at2(; ns = NS_MALLA_AT2)
    ps = [pico_at2(ParAT2(; MALLA_AT2..., N = n)) for n in ns]
    σ = [r.sigma for r in ps]
    difs = abs.(diff(σ))
    razones = difs[1:(end - 1)] ./ difs[2:end]
    (ns = collect(ns), picos = σ, difs = difs, razones = razones, ordenes = log2.(razones),
     orden = minimum(log2.(razones)), kkt = maximum(r.kkt for r in ps))
end

"""
    invariancia_escala_at2()

σ_pico·√ℓ con L/ℓ, h/ℓ y ancho/ℓ fijos: la ley ℓ^(-1/2) del AT2, medida sin ajustar ninguna recta.
`sin_raiz` es el control de mutación — sin el √ℓ los dos números difieren un 29 %.
"""
function invariancia_escala_at2()
    a = ParAT2(; L = 1.0, ℓ = 0.25, defecto = 0.4, ancho = 0.1, N = 80, fija_borde = true)
    b = ParAT2(; L = 2.0, ℓ = 0.50, defecto = 0.4, ancho = 0.2, N = 80, fija_borde = true)
    pa, pb = pico_at2(a).sigma, pico_at2(b).sigma
    (pico_a = pa, pico_b = pb, ell_a = a.ℓ, ell_b = b.ℓ,
     relativo = abs(pa * sqrt(a.ℓ) - pb * sqrt(b.ℓ)) / (pa * sqrt(a.ℓ)),
     sin_raiz = abs(pa - pb) / pa)
end

"""
    rama_homogenea_at2(q)

El control analítico del forward: en la barra homogénea con cotas naturales la rama uniforme es punto
estacionario exacto, así que el solver tiene que reproducir α(ε) = ε²/(ε²+ε_c²) y σ(ε) = Eε(1-α)², y
el pico σ_c = (3√3/16)√(EG_c/ℓ), a tolerancia de solver y no de discretización.
"""
function rama_homogenea_at2(q::ParAT2 = ParAT2(; L = 1.0, ℓ = 1.0, N = 40, defecto = 0.0))
    εs = RAMA_HOMOGENEA_AT2 .* eps_c(q)
    s = resolver_at2(εs .* q.L, q)
    eα = maximum(abs(maximum(s.alpha[k]) - alpha_homogenea(εs[k], q)) for k in eachindex(εs))
    eσ = maximum(abs(s.sigma[k] - sigma_homogenea(εs[k], q)) for k in eachindex(εs))
    pk = resolver_at2([eps_pico_homogenea(q) * q.L], q)
    # Hasta dónde vale la comparación, MEDIDO y no supuesto: la rama homogénea es un punto
    # estacionario a toda ε, pero deja de ser el mínimo que el esquema sigue. Con una rampa larga la
    # barra SALTA a daño total (α en amax) y desde ahí σ ya no es σ(ε) de la rama; medido el 14-nov
    # con ℓ = L: a 1,4 ε_c la rama se reproduce a 1e-13 y a 2,5 ε_c α vale 0,999999 contra 0,862 de la
    # rama. El criterio se mide en el rango pre-registrado de la prueba 4.1, y el salto se reporta.
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

"""
    gradiente_contra_fd_at2()

El gradiente adjunto contra diferencias finitas centradas en tres direcciones, donde el conjunto
activo no cambia: ÉSTA es la medición que decide si el objetivo 5 existe. `mutacion` es el gradiente
sin el término de rango uno del Hessiano reducido —el que saldría de «derivar el solver»—, que sale
con el signo cambiado.
"""
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

"""
    quiebre_at2()

El resultado del objetivo 5: en una etapa que SOSTIENE la carga anterior, α queda en la cota de
irreversibilidad con multiplicador cero, y las derivadas laterales del costo respecto de esa carga
difieren por un factor grande que NO se mueve al bajar τ cien veces. Es un quiebre, no truncamiento:
la aplicación carga→estado no es diferenciable, el adjunto da la derivada de la rama, y de ahí que en
el PMP clásico no esté disponible.
⚠️ CORREGIDO 14-nov: antes decía «en EDP solo se tenga C-estacionariedad». Atribución FALSA: ese
nombre es de Hintermüller-Kopacka 2009, no de Neitzel-Wick-Wollner. Ver docs/marco/objetivo5-phase-field.md § 4.2.
"""
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
    # Control: con la misma etapa 4 corrida un poco más arriba el conjunto activo queda vacío y el
    # mismo código concuerda con diferencias finitas. El desacuerdo es del punto, no del adjunto.
    limpia = [0.4, 0.6, 0.8, 0.85, 0.9] .* εc .* p.L
    r2 = gradiente_at2(limpia, p; objetivo = obj, omega = om)
    τ = 1e-6 * εc * p.L
    fd2 = (costo(limpia .+ τ .* e4) - costo(limpia .- τ .* e4)) / (2τ)
    # La separación del adjunto a CADA lateral, que es lo que NO depende del desempate de redondeo
    # del conjunto activo. Medido el 14-nov: `abs(adjunto - promedio) / promedio` pasa de 0,189 a
    # 0,049 al perturbar ℓ en 3,3e-12 relativo, o sea el criterio viejo —que pedía > 0,05— se pone
    # ROJO bajo una perturbación del tamaño del redondeo. Este no: su mínimo en las cinco corridas de
    # `docs/marco/re-medicion-impugnados/` es 1,507.
    sep = min(abs(r.grad[4] / taus[end]["derecha"]), abs(taus[end]["izquierda"] / r.grad[4]))
    (taus = taus, razones = razones, adjunto = r.grad[4], congelados = r.congelados,
     deriva_de_la_razon = maximum(razones) / minimum(razones) - 1,
     promedio = taus[end]["centrada"], separacion_laterales = sep,
     control_limpio = (congelados = r2.congelados, adjunto = r2.grad[4], fd = fd2,
                       error = abs(r2.grad[4] - fd2) / max(1.0, abs(fd2))))
end

"""
    seguimiento_at2(p, blancos; iter, v0)

El problema de control: llevar ∫α dx al programa pedido etapa por etapa. Con omega = 0 el óptimo es
cerrado —invertir Λ(d), que por rate-independencia es creciente—, así que el descenso proyectado sobre
los incrementos se mide contra un patrón conocido y no contra sí mismo.
"""
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
     # B2 del ataque a §9 (15-nov): `error_cerrado` es la precisión de la BISECCIÓN que construye la
     # referencia, no la distancia del descenso a ella. Esta es esa distancia, programa a programa.
     distancia_al_cerrado_pct = 100 * maximum(abs.(o.carga .- dex) ./ dex),
     J_cerrado = costo_at2(sx; objetivo = blancos, omega = 0.0).J,
     error_descenso = maximum(abs, o.sol.longitud .- blancos),
     relativo_descenso = maximum(abs, o.sol.longitud .- blancos) / blancos[end],
     pasos = length(o.historia), J0 = o.historia[1], J = o.J, historia = o.historia,
     costo_baja = all(diff(o.historia) .< 0), carga_monotona = all(diff(o.carga) .>= 0),
     congelados = g.congelados, fijos = count(nodos_fijos(p)), kkt = sx.kkt, p = p)
end

"""
    deficit_terminal_at2(p, blancos)

El descenso con blanco TERMINAL (escalar) contra el mismo blanco final del programa de seguimiento, a
omega = 0. Existe porque la afirmación «con blanco terminal el descenso se estanca» necesitaba una
corrida en el json y no la tenía: su única fuente era el docstring de `pesos_objetivo`.
⚠️ Medido el 14-nov ~21:50: NO se estanca. A omega = 0 el descenso ALCANZA el blanco terminal
(déficit 0,00 %) en 12-16 pasos, en las dos barras. Lo que el blanco terminal deja indeterminado es el
PROGRAMA, no el punto final — que es el argumento que de verdad sostiene el costo de seguimiento. El
estancamiento aparece con omega > 0, y entonces es del término de trabajo, no del blanco terminal.
"""
function deficit_terminal_at2(p::ParAT2, blancos::AbstractVector)
    K = length(blancos)
    term = blancos[end]
    dex = carga_exacta_at2(blancos, p)
    o = optimizar_at2(fill(dex[end] / K * 0.6, K), p; objetivo = term, omega = 0.0, iter = ITER_AT2)
    lam = o.sol.longitud[end]
    (blanco = term, longitud = lam, deficit_pct = 100 * (1 - lam / term),
     pasos = length(o.historia))
end

"""
    kkt_por_etapas_at2(; ks = KS_ETAPAS_AT2)

El residuo KKT de llegar a la MISMA carga final (0,9 ε_c L) en la MISMA barra `P_AT2_CONTROL`,
repartida en `K` etapas iguales. Existe por B3 del ataque a §9 (15-nov): la cifra que el informe
comparaba con el salto de una etapa era `orden_malla_at2().kkt`, que es el máximo sobre cuatro mallas
de una rampa de 30 etapas hasta 1,3 ε_c — otra barra y otra carga final. Medido acá, la frontera no
está en «una etapa contra varias»: con dos y con tres etapas el residuo sigue sobre 1e-2.
"""
function kkt_por_etapas_at2(; ks = KS_ETAPAS_AT2)
    p = P_AT2_CONTROL
    dfin = 0.9 * eps_c(p) * p.L
    Dict(string(K) => resolver_at2(collect(range(dfin / K, dfin; length = K)), p).kkt for K in ks)
end

"""
    objetivo5_at2()

El objetivo 5 completo como resultado publicable: los criterios del plan, la lectura del resultado, y
las series que dibuja `escribir_figura_at2_control`. Devuelve `(o, d)`: `o` va al json, `d` a la
figura, igual que `acuerdo()`.
"""
function objetivo5_at2()
    hom = rama_homogenea_at2()
    # M4 del ataque a §9 (15-nov): el default es ℓ_r = L, y el punto del salto depende de ℓ_r/L casi
    # por un factor 2. Esta segunda barra es la de la figura y la del experimento de malla.
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
        # Forward (tarea 4.1): el control analítico y la convergencia, con el criterio del plan.
        "la_rama_homogenea_reproduce_la_cerrada_bajo_1e9" =>
            max(hom.error_alpha, hom.error_sigma, hom.error_pico) < 1e-9,
        "el_orden_de_convergencia_en_malla_es_al_menos_085" => ml.orden >= 0.85,
        "la_invariancia_de_escala_pasa_1e9" => es.relativo < 1e-9,
        "mutacion_sin_el_raiz_de_ell_los_picos_difieren_mas_de_20pct" => es.sin_raiz > 0.2,
        # Adjunto (tarea 4.2): la medición que decide si el objetivo 5 existe, y su mutación.
        "el_gradiente_adjunto_concuerda_con_diferencias_finitas_bajo_1e5" =>
            maximum(gr.errores) < 1e-5,
        "mutacion_sin_el_termino_de_rango_uno_el_gradiente_cambia_de_signo" =>
            !gr.mutacion.mismo_signo,
        # El quiebre: el resultado.
        "las_derivadas_laterales_difieren_mas_de_3x" => minimum(qb.razones) > 3,
        "la_razon_de_las_laterales_no_se_mueve_1pct_al_bajar_tau_100x" =>
            qb.deriva_de_la_razon < 0.01,
        # ⚠️ CORREGIDO 14-nov ~21:30. Antes era
        #     "el_adjunto_no_es_el_promedio_de_las_laterales" =>
        #         abs(qb.adjunto - qb.promedio) > 0.05 * abs(qb.promedio)
        # y se ponía ROJO (0,0492 < 0,05) al perturbar ℓ en 3,3e-12 relativo: pinchaba el desempate de
        # redondeo del conjunto activo, no el quiebre. El de abajo exige que el adjunto esté separado
        # de LAS DOS laterales por al menos un 25 %, y eso aguanta las cinco corridas (mínimo 1,507).
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
        # Inestable (se mueve 29 % al perturbar ℓ en 3,3e-12 relativo): se registra, no se cita.
        "separacion_del_adjunto_a_la_lateral_mas_cercana" => qb.separacion_laterales,
        "nodos_congelados_en_la_etapa_que_sostiene" => qb.congelados[4],
        "nodos_de_la_barra_del_quiebre" => P_QUIEBRE_AT2.N + 1,
        # El blanco TERMINAL, que hasta el 14-nov ~21:50 el documento afirmaba «medido» sin corrida.
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
        # Límite medido, no criterio: la alternada necesita que la carga llegue en pasos chicos. Todo
        # este bloque es la MISMA barra y la MISMA carga final, repartida en K etapas (B3 del ataque
        # a §9): lo que manda es el tamaño del paso, no que haya más de una etapa.
        "kkt_si_una_sola_etapa_salta_a_09_epsc" => kkt[string(1)],
        "kkt_misma_barra_dos_etapas" => kkt[string(2)],
        "kkt_misma_barra_tres_etapas" => kkt[string(3)],
        "kkt_misma_barra_cinco_etapas" => kkt[string(5)],
        "kkt_misma_barra_diez_etapas" => kkt[string(10)],
        # ⚠️ NO es «la misma carga por etapas»: es el máximo sobre las cuatro mallas de la rampa de
        # 30 etapas hasta 1,3 ε_c. Otra barra, otra carga final. Así se llamaba y así engañaba.
        "kkt_de_la_rampa_de_30_etapas_peor_malla" => ml.kkt,
        "distancia_del_descenso_al_programa_cerrado_pct_con_defecto" => ct.distancia_al_cerrado_pct,
        "distancia_del_descenso_al_programa_cerrado_pct_homogenea" => hg.distancia_al_cerrado_pct,
        # Hasta dónde vale la rama homogénea como solución, no como punto estacionario.
        "la_rama_homogenea_calza_hasta_epsc" => hom.ultima_que_calza,
        "salta_a_dano_total_en_epsc" => hom.salto)

    o = Dict(
        "malla" => Dict("N" => ml.ns, "picos" => ml.picos, "diferencias" => ml.difs,
                        "razones" => ml.razones, "ordenes" => ml.ordenes,
                        # M10 del ataque a §9: §9.3 citaba `escala/ell_a` para describir ESTA barra.
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
