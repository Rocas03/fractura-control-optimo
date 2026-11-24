# Objetivo 2: el sistema linealizado en torno a una trayectoria nominal, que cuando importa es el óptimo del método
# directo. Aquí están la matriz fundamental Φ, los gramianos de controlabilidad y de observabilidad (midiendo solo la
# presión), los tests de Kalman en un instante, el control de energía mínima y un LQ que sigue al óptimo.
# Los gramianos se calculan por cuadratura sobre Φ y no con la EDO de Lyapunov que ya integra semilla.jl, para que la
# coincidencia con la semilla sea un cruce entre dos caminos distintos.

# ∂f/∂U: el caudal entra solo en la ecuación de P, dividido por el almacenamiento D.
∂fu(x, p) = SVector(0.0, 1 / D(x, p))

# C = [0 1]: en el pozo se mide la presión. El largo no se observa directamente, y eso es lo que cuantifica el
# gramiano de observabilidad.
const C_PRESION = SMatrix{1, 2}(0.0, 1.0)
const I2 = SMatrix{2, 2}(1.0, 0.0, 0.0, 1.0)

# Trayectoria nominal: el estado t -> x y el control t -> u.
struct Nominal{X, U}
    p::Par
    x::X
    u::U
    T::Float64
end

(nom::Nominal)(t) = nom.x(t)

# Nominal de un control dado como función de t, como las tres trayectorias de referencia de la semilla.
function nominal_control(p::Par, ufun, x0; T = p.T)
    s = integrar((y, t) -> f(y, ufun(t), p), SVector{2}(x0), (0.0, T))
    Nominal(p, τ -> SVector{2}(s(τ)), ufun, T)
end

# Nominal del óptimo del método directo: control constante a trozos y estado con la salida densa de los mismos
# tramos RK4 que usan avance y adjunto.
function nominal_optima(p::Par, q::AbstractVector)
    t = malla(p, length(q))
    tramos = avance(q, p)[2]
    idx(τ) = clamp(searchsortedlast(t, τ), 1, length(q))
    Nominal(p, τ -> x2(tramos[idx(τ)](clamp(τ, 0.0, p.T))), τ -> q[idx(τ)], p.T)
end

# A(t) = ∂f/∂x y B(t) = ∂f/∂U sobre la nominal.
matrices_lineales(nom::Nominal, t) = (x = nom.x(t); (∂f(x, nom.u(t), nom.p), ∂fu(x, nom.p)))

# Φ(·,t0) con salida densa: Φ' = A(t)Φ. Sirve hacia adelante y hacia atrás según el orden de (t0,t1).
function fundamental_densa(nom::Nominal, t0, t1)
    s = integrar((Y, t) -> vec(matrices_lineales(nom, t)[1] * SMatrix{2, 2}(Y)), vec(I2), (t0, t1))
    τ -> SMatrix{2, 2}(s(τ))
end

fundamental(nom::Nominal, t0, t1) = fundamental_densa(nom, t0, t1)(t1)

# Ψ(s) = Φ(T,s), de dΨ/ds = -Ψ A(s) con Ψ(T) = I, integrada de T a 0. Invertir Φ(s,0) daría lo mismo en teoría,
# pero arrastraría su condicionamiento.
function propagador_final(nom::Nominal, T = nom.T)
    s = integrar((Y, t) -> vec(-SMatrix{2, 2}(Y) * matrices_lineales(nom, t)[1]), vec(I2), (T, 0.0))
    τ -> SMatrix{2, 2}(s(τ))
end

simetrizar(W) = SMatrix{2, 2}((W + W') / 2)

# Cruce de los gramianos sobre la nominal del óptimo por las EDO de Lyapunov: W' = AW + WAᵀ + BBᵀ con W(0) = 0 y
# -V' = AᵀV + VA + CᵀC con V(T) = 0 dan W(T) y V(0) = M sin propagador ni cuadratura.
function cruce_lyapunov_nominal(nom::Nominal; T = nom.T)
    fW = function (Y, t)
        AB = matrices_lineales(nom, t)
        Wt = SMatrix{2, 2}(Y)
        vec(AB[1] * Wt + Wt * AB[1]' + AB[2] * AB[2]')
    end
    WL = simetrizar(SMatrix{2, 2}(integrar(fW, vec(0.0 * I2), (0.0, T))(T)))
    fV = function (Y, t)
        AB = matrices_lineales(nom, t)
        Vt = SMatrix{2, 2}(Y)
        vec(-(AB[1]' * Vt + Vt * AB[1] + C_PRESION' * C_PRESION))
    end
    ML = simetrizar(SMatrix{2, 2}(integrar(fV, vec(0.0 * I2), (T, 0.0))(0.0)))
    λWq = eigvals(Symmetric(Matrix(gramiano_controlabilidad(nom; T))))
    λMq = eigvals(Symmetric(Matrix(gramiano_observabilidad(nom; T))))
    λWl = eigvals(Symmetric(Matrix(WL)))
    λMl = eigvals(Symmetric(Matrix(ML)))
    eW = maximum(abs.(λWl - λWq)) / max(1.0, maximum(abs, λWq))
    eM = maximum(abs.(λMl - λMq)) / max(1.0, maximum(abs, λMq))
    (lambda_W_lyapunov = λWl, lambda_M_lyapunov = λMl,
     error_rel_controlabilidad = eW, error_rel_observabilidad = eM, error_rel_maximo = max(eW, eM))
end

# W = ∫ Φ(T,s) B Bᵀ Φ(T,s)ᵀ ds y M = ∫ Φ(s,0)ᵀ Cᵀ C Φ(s,0) ds, por cuadratura de Gauss-Kronrod.
function gramiano_controlabilidad(nom::Nominal; T = nom.T)
    Ψ = propagador_final(nom, T)
    W = quadgk(s -> (g = Ψ(s) * ∂fu(nom.x(s), nom.p); g * g'), 0.0, T; rtol = 1e-12, atol = 1e-14)[1]
    SMatrix{2, 2}(W)
end

function gramiano_observabilidad(nom::Nominal; T = nom.T)
    Φ = fundamental_densa(nom, 0.0, T)
    M = quadgk(s -> (g = C_PRESION * Φ(s); g' * g), 0.0, T; rtol = 1e-12, atol = 1e-14)[1]
    SMatrix{2, 2}(M)
end

# Rango numérico de una 2×2 por valores singulares, con umbral relativo √eps.
function _rango_sv(K)
    sv = svdvals(K)
    (rango = count(>(sqrt(eps(Float64)) * max(sv[1], 1.0)), sv), sv = SVector{2}(sv))
end

# Rango de [B AB] y de [C; CA] con las matrices congeladas en un instante t.
kalman_congelado(nom::Nominal, t) = (AB = matrices_lineales(nom, t); _rango_sv(hcat(AB[2], AB[1] * AB[2])))
observabilidad_congelada(nom::Nominal, t) =
    (A = matrices_lineales(nom, t)[1]; _rango_sv(vcat(C_PRESION, C_PRESION * A)))

# Control de energía mínima que lleva z(0) = 0 a z(T) = z_T en el sistema linealizado:
# u(s) = B(s)'Ψ(s)' W⁻¹ z_T, y ∫₀ᵀ u² = z_T' W⁻¹ z_T. Falla, como debe, si W es singular.
function energia_minima(nom::Nominal, zT; T = nom.T)
    W = gramiano_controlabilidad(nom; T)
    z = SVector{2}(zT)
    c = SVector{2}(Symmetric(Matrix(W)) \ Vector(z))
    Ψ = propagador_final(nom, T)
    (u = s -> dot(∂fu(nom.x(s), nom.p), Ψ(s)' * c), energia = dot(z, c))
end

# LQ de seguimiento: -Ṗ = PA + A'P - PB R⁻¹ B'P + Q, P(T) = Qf, integrada hacia atrás con salida densa.
# La ganancia es K(t) = -R⁻¹B(t)'P(t), y el control de seguimiento q(t) = q̄(t) + K(t)(x(t) - x̄(t)).
function riccati(nom::Nominal; Q = I2, R = 1.0, Qf = I2, T = nom.T)
    s = integrar(vec(SMatrix{2, 2}(Qf)), (T, 0.0)) do Y, t
        P = simetrizar(SMatrix{2, 2}(Y))
        A, B = matrices_lineales(nom, t)
        vec(-(P * A + A' * P - (P * B) * (B' * P) / R + SMatrix{2, 2}(Q)))
    end
    τ -> simetrizar(SMatrix{2, 2}(s(τ)))
end

ganancia_lq(nom::Nominal, P, t; R = 1.0) = -(matrices_lineales(nom, t)[2]' * P(t)) / R

# Simula el sistema no lineal con el control de seguimiento, recortado a [0, qmax] porque es un caudal. Con x0 o p
# distintos de los nominales es el experimento de perturbación; con lazo_abierto se aplica solo el control nominal.
function lazo_cerrado(nom::Nominal, P; x0 = nom.x(0.0), p = nom.p, R = 1.0, T = nom.T, lazo_abierto = false)
    control(t, x) = begin
        q̄ = nom.u(t)
        lazo_abierto ? clamp(q̄, 0.0, p.qmax) :
            clamp(q̄ + dot(ganancia_lq(nom, P, t; R), SVector{2}(x) - nom.x(t)), 0.0, p.qmax)
    end
    s = integrar((y, t) -> (x = x2(y); push(f(x, control(t, x), p), L(x, control(t, x), p))),
                 SVector(x0[1], x0[2], 0.0), (0.0, T))
    (sol = s, x = τ -> x2(s(τ)), u = τ -> control(τ, x2(s(τ))), xT = x2(s.u[end]), costo = s.u[end][3],
     J = Φ(x2(s.u[end]), p) + s.u[end][3])
end

# Resultado del objetivo 2, que regenerar_final.jl escribe en final/resultados/objetivo2_lineal.json.

const N_OBJ2 = 50                                  # tramos del control del óptimo nominal

# Los dos casos del experimento LQ (P0 o η un 20 % mayores), en un solo lugar para que el json y la figura corran
# exactamente lo mismo.
casos_lq(p) = (("P0", SVector(p.ℓ0, (1 + PERT) * p.P0), p),
               ("eta", SVector(p.ℓ0, p.P0), con(p; η = (1 + PERT) * p.η)))

# Series temporales del experimento LQ para la figura, que el json no guarda: por caso, la nominal, el lazo abierto
# perturbado y el LQ perturbado, con sus controles.
function datos_figura_lq(p::Par = Par(); N = N_OBJ2, npts = 401)
    nom = nominal_optima(p, optimizar(p, N).x)
    Q = SMatrix{2, 2}(0.0, 0.0, 0.0, p.wp)
    Qf = SMatrix{2, 2}(p.wT, 0.0, 0.0, 0.0)
    P = riccati(nom; Q, R = p.wq, Qf)
    ts = collect(range(0.0, p.T, length = npts))
    casos = map(casos_lq(p)) do (cual, x0, pp)
        ab = lazo_cerrado(nom, P; x0, p = pp, R = p.wq, lazo_abierto = true)
        ce = lazo_cerrado(nom, P; x0, p = pp, R = p.wq)
        (cual = cual, ab = (ell = [ab.x(t)[1] for t in ts], u = [ab.u(t) for t in ts], xT = ab.xT),
         ce = (ell = [ce.x(t)[1] for t in ts], u = [ce.u(t) for t in ts], xT = ce.xT))
    end
    (t = ts, nominal = (ell = [nom.x(t)[1] for t in ts], u = [clamp(nom.u(t), 0.0, p.qmax) for t in ts]),
     casos = casos, ellT_nominal = nom.x(p.T)[1], elld = p.ℓd, qmax = p.qmax, pert = PERT)
end
const TS_BARRIDO = (0.5, 1.0, 1.5, 2.0)            # horizontes del barrido de condicionamiento
const ETAS_BARRIDO = (0.0, 0.2, 0.5, 1.0)          # filtración
const Z_BLANCO = SVector(0.013, -0.021)            # blanco del control de energía mínima
const PERT = 0.2                                   # perturbación relativa del experimento LQ
const NPTS_NOMINAL = 201                           # muestreo de la nominal para buscar la degeneración
const TOL_SEMILLA = 1e-8                           # criterio del contraste cruzado con la semilla
const UMBRAL_DEGEN = 1e-10                         # σ₂/σ₁ bajo el cual se declara degeneración

filas(A) = [collect(r) for r in eachrow(A)]
cond2(λ) = (m = minimum(λ); m ≤ 0 ? Inf : maximum(λ) / m)

# Las tres nominales de los gramianos de referencia de semilla.jl: la presurizada, con u(t) = 0.5 + 0.1 sin t desde
# (1, 0.3), y las dos de reposo.
nominales_semilla(p = Par()) = Dict(
    "trayectoria_presurizada" => nominal_control(con(p; n = 2.0), t -> 0.5 + 0.1 * sin(t), SVector(1.0, 0.3); T = 2.0),
    "reposo_n2" => nominal_control(con(p; n = 2.0), _ -> 0.0, SVector(1.0, 0.0); T = 2.0),
    "reposo_n1" => nominal_control(con(p; n = 1.0), _ -> 0.0, SVector(1.0, 0.0); T = 2.0))

# Gramianos por cuadratura contra los de la semilla, en esas tres nominales.
function contraste_semilla(sem = semilla_lectura(), p = Par())
    out = Dict{String, Any}()
    peor = 0.0
    for (caso, nom) in nominales_semilla(p)
        ref = sem["gramianos_$caso"]
        mW = eigvals(Symmetric(Matrix(gramiano_controlabilidad(nom; T = 2.0))))
        mM = eigvals(Symmetric(Matrix(gramiano_observabilidad(nom; T = 2.0))))
        eW = maximum(abs.(mW - ref["controlabilidad"])) / max(1.0, maximum(abs, ref["controlabilidad"]))
        eM = maximum(abs.(mM - ref["observabilidad"])) / max(1.0, maximum(abs, ref["observabilidad"]))
        peor = max(peor, eW, eM)
        out[caso] = Dict{String, Any}("controlabilidad_cuadratura" => mW, "controlabilidad_semilla" => ref["controlabilidad"],
                                      "observabilidad_cuadratura" => mM, "observabilidad_semilla" => ref["observabilidad"],
                                      "error_rel_controlabilidad" => eW, "error_rel_observabilidad" => eM,
                                      # El n de este caso, no el de Par().
                                      "n" => nom.p.n)
    end
    out["error_rel_maximo"] = peor
    out["metodo"] = "la semilla integra la EDO de Lyapunov W' = AW + WA' + BB' y ∫X'C'CX junto al estado; " *
                    "acá es cuadratura de Gauss–Kronrod sobre el propagador Ψ(s) = Φ(T,s) y Φ(s,0). " *
                    "Son dos caminos distintos, no la misma cuenta dos veces."
    out
end

# Nominal óptima, gramianos y su condicionamiento, degeneración en P = 0, energía mínima y LQ de seguimiento.
function objetivo2(p::Par = Par(); N = N_OBJ2)
    opt = optimizar(p, N)
    nom = nominal_optima(p, opt.x)
    W, M = gramiano_controlabilidad(nom), gramiano_observabilidad(nom)
    λW = eigvals(Symmetric(Matrix(W)))
    λM = eigvals(Symmetric(Matrix(M)))

    # Lo interesante no es que sea controlable, sino cuán mal condicionado está: eso dice qué dirección
    # del estado cuesta caro mover, y cuánto ruido de presión se amplifica al reconstruir ℓ.
    barrido = Dict{String, Any}[]
    for T in TS_BARRIDO, η in ETAS_BARRIDO
        pb = con(p; η)
        nb = nominal_optima(pb, optimizar(pb, N).x)
        lw = eigvals(Symmetric(Matrix(gramiano_controlabilidad(nb; T))))
        lm = eigvals(Symmetric(Matrix(gramiano_observabilidad(nb; T))))
        push!(barrido, Dict{String, Any}("T" => T, "eta" => η, "lambda_W" => lw, "lambda_M" => lm,
                                         "cond_W" => cond2(lw), "cond_M" => cond2(lm)))
    end

    # La degeneración que el modelo predice en P = 0: con n = 2, ∂F/∂ℓ = ∂F/∂P = 0 ahí, la primera fila
    # de A se anula y ℓ queda desacoplada. Se mide sobre la nominal del óptimo (donde P > 0) y en reposo.
    ts = range(0.0, p.T, length = NPTS_NOMINAL)
    rk = [kalman_congelado(nom, t) for t in ts]
    ro = [observabilidad_congelada(nom, t) for t in ts]
    razon_k = [r.sv[2] / r.sv[1] for r in rk]
    razon_o = [r.sv[2] / r.sv[1] for r in ro]
    Pmin = minimum(nom.x(t)[2] for t in ts)
    nsem = nominales_semilla(p)
    kr = kalman_congelado(nsem["reposo_n2"], 1.0)
    orp = observabilidad_congelada(nsem["reposo_n2"], 1.0)
    degenerada_en_la_nominal = minimum(razon_k) < UMBRAL_DEGEN
    degen = Dict{String, Any}(
        "P_min_en_la_nominal" => Pmin,
        "n_reposo" => nsem["reposo_n2"].p.n,
        "kalman_razon_sv_min" => minimum(razon_k), "kalman_rango_min" => minimum(r.rango for r in rk),
        "observabilidad_razon_sv_min" => minimum(razon_o), "observabilidad_rango_min" => minimum(r.rango for r in ro),
        "reposo_kalman_rango" => kr.rango, "reposo_kalman_razon_sv" => kr.sv[2] / kr.sv[1],
        "reposo_observabilidad_rango" => orp.rango, "reposo_observabilidad_razon_sv" => orp.sv[2] / orp.sv[1],
        "degenerada_en_la_nominal" => degenerada_en_la_nominal,
        "lectura" => degenerada_en_la_nominal ?
            "La degeneración es real sobre la nominal: hay puntos con σ₂/σ₁ < $UMBRAL_DEGEN." :
            "Sobre la nominal del óptimo no hay degeneración (P > 0 en todo el horizonte): controlabilidad y " *
            "observabilidad uniformes, con cota σ₂/σ₁ ≥ $(minimum(razon_k)) y $(minimum(razon_o)). La degeneración " *
            "existe y vive en P = 0, donde el rango cae a 1 en los dos tests: en reposo ℓ no se controla ni se observa.")

    cruce = cruce_lyapunov_nominal(nom)
    em = energia_minima(nom, Z_BLANCO)
    zT = integrar((y, t) -> let (A, B) = matrices_lineales(nom, t); A * y + B * em.u(t) end,
                  SVector(0.0, 0.0), (0.0, p.T)).u[end]
    energia_cuadratura = quadgk(s -> em.u(s)^2, 0.0, p.T; rtol = 1e-12)[1]
    forma = dot(Z_BLANCO, Symmetric(Matrix(W)) \ Vector(Z_BLANCO))

    # LQ de seguimiento con los pesos que salen del problema: Qf = ∂²Φ = diag(w_T, 0), Q = ∂²L = diag(0, w_p),
    # R = w_q. Es Gauss–Newton del costo original en torno al óptimo, no pesos elegidos a gusto.
    Q = SMatrix{2, 2}(0.0, 0.0, 0.0, p.wp)
    Qf = SMatrix{2, 2}(p.wT, 0.0, 0.0, 0.0)
    P = riccati(nom; Q, R = p.wq, Qf)
    ℓT̄ = nom.x(p.T)[1]
    exp_lq = Dict{String, Any}()
    factores = Float64[]
    for (cual, x0, pp) in casos_lq(p)
        ab = lazo_cerrado(nom, P; x0, p = pp, R = p.wq, lazo_abierto = true)
        ce = lazo_cerrado(nom, P; x0, p = pp, R = p.wq)
        da, dc = abs(ab.xT[1] - ℓT̄), abs(ce.xT[1] - ℓT̄)
        push!(factores, da / dc)
        exp_lq[cual] = Dict{String, Any}(
            "perturbacion_relativa" => PERT,
            "desviacion_terminal_lazo_abierto" => da, "desviacion_terminal_LQ" => dc, "factor" => da / dc,
            "error_a_ell_d_lazo_abierto" => abs(ab.xT[1] - p.ℓd), "error_a_ell_d_LQ" => abs(ce.xT[1] - p.ℓd),
            "error_a_ell_d_nominal" => abs(ℓT̄ - p.ℓd),
            "J_lazo_abierto" => ab.J, "J_LQ" => ce.J, "ell_T_lazo_abierto" => ab.xT[1], "ell_T_LQ" => ce.xT[1],
            "U_LQ_min" => minimum(ce.u(t) for t in range(0.0, p.T, length = 401)),
            "U_LQ_max" => maximum(ce.u(t) for t in range(0.0, p.T, length = 401)),
            "recorte_activo" => any(t -> (v = ce.u(t); v <= 1e-12 || v >= pp.qmax - 1e-12),
                                    range(0.0, p.T, length = 401)))
    end
    nulo = lazo_cerrado(nom, P; R = p.wq)

    # Con Qf = diag(wT, 0), P(T) es singular por construcción y P no puede ser definida positiva en todo [0, T].
    # Se verifica P ⪰ 0, la condición de borde exacta y el residuo de la propia EDO.
    tsP = range(0.0, p.T, length = 21)
    λminP = [minimum(eigvals(Symmetric(Matrix(P(t))))) for t in tsP]
    borde_P = maximum(abs.(P(p.T) - Qf))
    hR = 1e-5
    residual_P = maximum(range(hR, p.T - hR, length = 50)) do t
        A, B = matrices_lineales(nom, t)
        Pt = P(t)
        maximum(abs.((P(t + hR) - P(t - hR)) / 2hR + (Pt * A + A' * Pt - (Pt * B) * (B' * Pt) / p.wq + Q)))
    end

    criterios = Dict{String, Any}(
        "gramianos_contra_la_semilla" => contraste_semilla(semilla_lectura(), p)["error_rel_maximo"] < TOL_SEMILLA,
        "los_gramianos_de_la_nominal_coinciden_por_lyapunov" => cruce.error_rel_maximo < TOL_SEMILLA,
        "W_simetrica_definida_positiva" => maximum(abs.(W - W')) < 1e-12 * max(1.0, maximum(abs, W)) && minimum(λW) > 0,
        "M_simetrica_definida_positiva" => maximum(abs.(M - M')) < 1e-12 * max(1.0, maximum(abs, M)) && minimum(λM) > 0,
        "energia_minima_alcanza_el_blanco" => maximum(abs.(zT - Z_BLANCO)) < 1e-8 * maximum(abs, Z_BLANCO),
        # Esto no es un cruce: energia_minima devuelve la misma expresión zᵀW⁻¹z y coincide bit a bit. El cruce con
        # contenido es el siguiente, contra la cuadratura del control.
        "energia_y_forma_cuadratica_identidad_no_cruce" => abs(em.energia - forma) < 1e-9 * abs(em.energia),
        "energia_igual_a_la_cuadratura_del_control" => abs(energia_cuadratura - em.energia) < 1e-8 * abs(em.energia),
        "degeneracion_en_reposo_rango_1" => kr.rango == 1 && orp.rango == 1,
        "riccati_semidefinida_positiva" => minimum(λminP) > -1e-10 * maximum(abs, P(0.0)),
        "riccati_definida_positiva_antes_de_T" => all(λminP[1:end-1] .> 0),
        "riccati_condicion_de_borde_exacta" => borde_P == 0.0,
        "riccati_residual_de_su_edo" => residual_P < 1e-8,
        "perturbacion_nula_reproduce_la_nominal" => maximum(abs.(nulo.xT - nom.x(p.T))) < 1e-9,
        "LQ_reduce_la_desviacion_factor_3" => all(f -> f > 3, factores))

    Dict{String, Any}(
        "parametros" => comodict(p), "N_tramos" => N, "nominal" => Dict{String, Any}(
            "J" => opt.fun, "exito_ipopt" => opt.exito, "status" => string(opt.status), "iteraciones" => opt.iter,
            "cotas_fijadas" => opt.fijadas, "q" => opt.x, "ell_T" => nom.x(p.T)[1], "P_T" => nom.x(p.T)[2]),
        "gramiano_controlabilidad" => Dict{String, Any}("W" => filas(W), "lambda" => λW, "cond" => cond2(λW)),
        # Gramianos por Lyapunov sobre la nominal del óptimo, que es la que se reporta.
        "cruce_lyapunov_en_la_nominal" => Dict{String, Any}(
            "lambda_W_lyapunov" => cruce.lambda_W_lyapunov,
            "lambda_M_lyapunov" => cruce.lambda_M_lyapunov,
            "error_rel_controlabilidad" => cruce.error_rel_controlabilidad,
            "error_rel_observabilidad" => cruce.error_rel_observabilidad,
            "error_rel_maximo" => cruce.error_rel_maximo,
            "metodo" => "W' = AW + WA' + BB' con W(0) = 0 y -V' = A'V + VA + C'C con V(T) = 0, " *
                        "integradas sobre la nominal del óptimo; sin propagador y sin cuadratura. " *
                        "El acuerdo es tres órdenes peor que en las tres trayectorias suaves de la " *
                        "semilla. Atribución medida aparte, que este script NO calcula (está en " *
                        "docs/marco/ataque-s5/CRUCE-LYAPUNOV-MEDIDO.md): el origen es que A(τ) es " *
                        "discontinua en los nudos del control, pero el mecanismo no es que un paso " *
                        "pise un salto —poner tstops en los 49 nudos recupera 6,1x de los 644x, y los " *
                        "mismos 49 tstops sobre una nominal suave compran 9,6x por sí solos—, sino el " *
                        "número de pasos (1533 contra 27 en la suave), y es artefacto de tolerancia " *
                        "(apretar reltol de 1e-12 a 1e-14 baja la cifra 57x). Escala e interpolante " *
                        "del estado quedan refutados. El mecanismo se apoya en dos indicios indirectos, " *
                        "no en una cuenta del error local por paso.",
            "tolerancia" => TOL_SEMILLA,
            "procedencia_de_la_tolerancia" => "heredada de TOL_SEMILLA, fijada antes para el " *
                                              "contraste con la semilla; no se eligió viendo este número"),
        "gramiano_observabilidad" => Dict{String, Any}("M" => filas(M), "lambda" => λM, "cond" => cond2(λM),
                                                       "C" => [0.0, 1.0], "lectura" => "C = [0 1]: se mide solo la presión en el pozo."),
        "barrido_condicionamiento" => barrido, "degeneracion" => degen,
        "energia_minima" => Dict{String, Any}("z_T" => Z_BLANCO, "energia" => em.energia,
                                              "forma_cuadratica_zWz" => forma, "energia_por_cuadratura" => energia_cuadratura,
                                              "z_T_alcanzado" => zT),
        "lq_seguimiento" => Dict{String, Any}("Q" => filas(Q), "Qf" => filas(Qf), "R" => p.wq,
                                              "pesos" => "segunda variación del costo original (Gauss–Newton): Qf = ∂²Φ, Q = ∂²L, R = w_q",
                                              "P_0" => filas(P(0.0)), "lambda_P_0" => eigvals(Symmetric(Matrix(P(0.0)))),
                                              "lambda_min_en_la_grilla" => λminP, "borde_max_abs_P_T_menos_Qf" => borde_P,
                                              "residual_edo_max" => residual_P,
                                              "experimento" => exp_lq),
        # P definida positiva en todo [0, T] no se cumple, por la estructura del costo: queda registrado con su razón
        # y lo reemplazan los criterios de arriba.
        "criterios_retirados" => Dict{String, Any}("riccati_definida_positiva" => Dict{String, Any}(
            "cumplido" => all(λminP .> 0), "fijado_en" => "plan 2026-11-14, F1 tarea 1.4 paso 1",
            "razon" => "Imposible con los pesos del problema, y por estructura del costo, no por el código: " *
                       "Qf = ∂²Φ = diag(w_T, 0) no penaliza la presión, así que P(T) = Qf es singular por diseño. " *
                       "Medido: λ_min(P(t)) > 0 para t < T, decayendo linealmente (9,99e-4 en t = 1,999), " *
                       "λ_min(P(T)) = 0 y max|P(T) − Qf| = 0,0 exacto. Con Q = Qf = I, donde la pregunta sí tiene " *
                       "sentido, P sí es definida positiva: lo prueba el testset «1.4 Riccati» de test/pruebas.jl.",
            "reemplazado_por" => ["riccati_semidefinida_positiva", "riccati_definida_positiva_antes_de_T",
                                  "riccati_condicion_de_borde_exacta", "riccati_residual_de_su_edo"])),
        "contraste_semilla" => contraste_semilla(semilla_lectura(), p),
        "criterios" => criterios)
end
