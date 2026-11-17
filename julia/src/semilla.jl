# Port de semilla/verificar_modelo.py (e31f323), sección por sección y con los mismos assert. No resuelve el proyecto.
# Las 8 identidades se comprueban en números, no en símbolos: 20 puntos al azar en BigFloat de 256 bits, derivadas por
# paso complejo. Lo demás usa f y ∂f de «Parámetros y modelo», que son el rhs y el jac de la semilla de Python.

const PASO_COMPLEJO_BIG = big"1e-60"
dcomplejo(g, x) = imag(g(complex(x, PASO_COMPLEJO_BIG))) / PASO_COMPLEJO_BIG

const IDENTIDADES = ("balance_almacenamiento", "jacobiano_F_L", "jacobiano_F_p", "jacobiano_h_L", "jacobiano_h_p",
                     "kalman_congelado", "PMP_estacionariedad", "HJB_minimizacion")

# Devuelve, por identidad, el máximo de |expresión| en los puntos.
function identidades_bigfloat(rng; puntos = 20)
    peor = Dict(k => big(0.0) for k in IDENTIDADES)
    setprecision(BigFloat, 256) do
        positivo() = big"0.3" + big"1.2" * rand(rng, BigFloat)
        signo() = 2 * rand(rng, BigFloat) - 1
        for k in 1:puntos
            n = BigFloat((1, 2, 3, 5)[mod1(k, 4)])
            Lv, pv, qv, a, b, C, λ, wp, wq = (positivo() for _ in 1:9)
            z1, z2, vL, vp = (signo() for _ in 1:4)
            Dd(L) = C + b * L^2
            Fd(L, p) = a * L^(n / 2) * p^n
            Rd(L, p) = 2λ * L * p + 2b * L * p * Fd(L, p)
            hd(L, p, q) = (q - Rd(L, p)) / Dd(L)
            D_, F, R, h = Dd(Lv), Fd(Lv, pv), Rd(Lv, pv), hd(Lv, pv, qv)
            FL, Fp = n * F / (2Lv), n * F / pv
            RL, Rp = 2λ * pv + b * (n + 2) * pv * F, 2λ * Lv + 2b * (n + 1) * Lv * F
            hL, hp = -RL / D_ - 2b * Lv * h / D_, -Rp / D_
            A = [FL Fp; hL hp]
            B = [zero(D_), 1 / D_]
            Hq(q) = z1 * F + z2 * hd(Lv, pv, q) - (wp * pv^2 + wq * q^2) / 2
            qstar = z2 / (wq * D_)
            qhjb = -vp / (wq * D_)
            Hb = F * vL - R * vp / D_ + wp * pv^2 / 2 - vp^2 / (2wq * D_^2)
            valores = (
                dcomplejo(L -> Dd(L) * pv, Lv) * F + dcomplejo(p -> Dd(Lv) * p, pv) * h - (qv - 2λ * Lv * pv),
                dcomplejo(L -> Fd(L, pv), Lv) - FL,
                dcomplejo(p -> Fd(Lv, p), pv) - Fp,
                dcomplejo(L -> hd(L, pv, qv), Lv) - hL,
                dcomplejo(p -> hd(Lv, p, qv), pv) - hp,
                det([B A * B]) + Fp / D_^2,
                dcomplejo(Hq, qstar),
                (wp * pv^2 + wq * qhjb^2) / 2 + F * vL + hd(Lv, pv, qhjb) * vp - Hb)
            for (nombre, v) in zip(IDENTIDADES, valores)
                peor[nombre] = max(peor[nombre], abs(v))
            end
        end
    end
    peor
end

# NoSpecialize: Vern9 se compila una vez por tipo de estado y no una por campo (cada clausura es un tipo nuevo); medido
# en T6: 7,4 s de compilación en frío con la especialización por defecto.
const SIN_ESPECIALIZAR = parentmodule(ODEProblem).NoSpecialize
function integrar_semilla(g, y0, tspan)
    s = solve(ODEProblem{false, SIN_ESPECIALIZAR}((y, _, t) -> g(y, t), y0, tspan), Vern9(); reltol = 2e-11,
              abstol = 2e-13)
    s.retcode == ReturnCode.Success || error("integrar_semilla: $(s.retcode)")
    s
end

function semilla_lectura()
    res = Dict{String, Any}()
    res["versions"] = merge(Dict("julia" => string(VERSION)),
                            Dict(string(m) => string(pkgversion(m)) for m in (OrdinaryDiffEqVerner, NLPModels,
                                                                               NLPModelsIpopt, QuadGK, JSON, StaticArrays)))

    # Identidades dimensionales, con L, p > 0 para derivar potencias reales.
    peor = identidades_bigfloat(Xoshiro(4703))
    for k in IDENTIDADES
        @printf("   %-24s max |expr| = %.1e\n", k, Float64(peor[k]))
    end
    ids = Dict{String, Any}(k => peor[k] < big"1e-50" for k in IDENTIDADES)
    @assert all(ids[k] for k in IDENTIDADES)
    ids["metodo"] = "numérico: 20 puntos al azar, BigFloat de 256 bits, derivadas por paso complejo; " *
                    "no es una verificación simbólica"
    res["identidades_numericas_bigfloat"] = ids

    # Jacobiano analítico contra paso complejo, en Float64.
    rng = Xoshiro(4703)
    errores = Float64[]
    for n in (1.0, 2.0, 3.5), _ in 1:10
        x = SVector{2}(0.3 .+ 1.2 .* rand(rng, 2))
        u = 0.1 + 0.9 * rand(rng)
        p = con(Par(); n)
        numerico = hcat((imag.(f(x + 1e-25im * e, u, p)) ./ 1e-25 for e in (SVector(1.0, 0.0), SVector(0.0, 1.0)))...)
        push!(errores, maximum(abs, ∂f(x, u, p) - numerico))
    end
    res["jacobiano_paso_complejo_error_max"] = maximum(errores)
    @assert maximum(errores) < 1e-11

    # Coordenadas (ℓ, P) contra (ℓ, S, fuga) y balance de masa; n = 2, c = 1, η = 0,2 (Par por defecto).
    p = Par()
    T, x0 = 2.0, SVector(1.0, 0.3)
    control(t) = 0.5 + 0.1 * sin(t)
    Bu(x, q) = SVector(0.0, 1 / D(x, q))
    nom = integrar_semilla((x, t) -> f(x, control(t), p), x0, (0.0, T))
    function almacenamiento(y, t)
        ℓ, S = y[1], y[2]
        P = S / (1 + ℓ^2)
        SVector(ℓ * P^2, control(t) - 0.2ℓ * P, 0.2ℓ * P)
    end
    st = integrar_semilla(almacenamiento, SVector(1.0, 0.6, 0.0), (0.0, T))
    malla = range(0.0, T, length = 1001)
    ys = st.(malla)
    err_coord = maximum(maximum(abs, nom(t) - SVector(y[1], y[2] / (1 + y[1]^2))) for (t, y) in zip(malla, ys))
    err_masa = maximum(abs(y[2] + y[3] - 0.6 - (0.5t + 0.1 * (1 - cos(t)))) for (t, y) in zip(malla, ys))
    res["coordenadas_y_balance"] = Dict{String, Any}("error_max_estados" => err_coord, "error_max_balance" => err_masa,
                                                     "estado_final" => collect(nom.u[end]))
    @assert err_coord < 1e-8 && err_masa < 1e-8

    # Referencia independiente: q = η = 0, n = 1, S constante y solución implícita.
    pa = con(p; n = 1.0, η = 0.0)
    apagado = integrar_semilla((x, _) -> f(x, 0.0, pa), x0, (0.0, 3.0))
    ta = range(0.0, 3.0, length = 701)
    xa = apagado.(ta)
    ℓa, Pa, S0 = first.(xa), last.(xa), 0.6
    resid = @. 2 * (sqrt(ℓa) - 1) + 2 / 5 * (ℓa^2.5 - 1) - S0 * ta
    res["apagado_sin_filtracion"] = Dict{String, Any}("residuo_solucion_implicita" => maximum(abs, resid),
                                                      "error_S_constante" => maximum(abs, @. (1 + ℓa^2) * Pa - S0),
                                                      "L_final" => ℓa[end])
    @assert maximum(abs, resid) < 1e-8

    # Adimensionalización comprobada con escalas no unitarias.
    E, KI, L0, v0, cc, ll, nn = 8.0, 2.0, 1.3, 0.2, 0.7, 0.04, 2.5
    bb = 2π / E
    ps, ts = KI / sqrt(π * L0), L0 / v0
    qs = bb * L0^2 * ps / ts
    LL, pp, qq = L0 * 0.9, ps * 0.6, qs * 0.4
    vel = v0 * (pp * sqrt(π * LL) / KI)^nn
    fisico = SVector(vel, (qq - 2ll * LL * pp - 2bb * pp * LL * vel) / (cc + bb * LL^2))
    escalado = fisico .* ts ./ SVector(L0, ps)
    nd = f(SVector(0.9, 0.6), 0.4, con(p; n = nn, c = cc / (bb * L0^2), η = 2ll * ts / (bb * L0)))
    res["adimensionalizacion_error"] = maximum(abs, escalado - nd)
    @assert maximum(abs, escalado - nd) < 1e-12

    # La aproximación lineal debe tener residuo de orden ε².
    dir(t) = cos(1.5t)
    vari = integrar_semilla((z, t) -> (x = nom(t); ∂f(x, control(t), p) * z + Bu(x, p) * dir(t)),
                            SVector(0.0, 0.0), (0.0, T))
    lin = Dict{String, Any}[]
    for ε in (0.04, 0.02, 0.01, 0.005)
        pert = integrar_semilla((x, t) -> f(x, control(t) + ε * dir(t), p), x0, (0.0, T))
        e = norm(pert.u[end] - nom.u[end] - ε * vari.u[end])
        push!(lin, Dict{String, Any}("eps" => ε, "error" => e, "error_div_eps2" => e / ε^2))
    end
    res["linealizacion"] = lin
    @assert all(3.7 < lin[i]["error"] / lin[i+1]["error"] < 4.3 for i in 1:3)

    # Adjunto con H = ζ·f − costo: derivada del costo frente a diferencias finitas.
    wT, wP, wU, ℓdes = 10.0, 1.0, 0.2, 1.7
    function objetivo_semilla(ε)
        campo(y, t) = (u = control(t) + ε * dir(t); push(f(x2(y), u, p), (wP * y[2]^2 + wU * u^2) / 2))
        s = integrar_semilla(campo, SVector(x0[1], x0[2], 0.0), (0.0, T))
        wT * (s.u[end][1] - ℓdes)^2 / 2 + s.u[end][3]
    end
    function adjunto_semilla(y, t)
        x = nom(t)
        z = x2(y)
        dz = -∂f(x, control(t), p)' * z + SVector(0.0, wP * x[2])
        push(dz, (wU * control(t) - z ⋅ Bu(x, p)) * dir(t))
    end
    adj = integrar_semilla(adjunto_semilla, SVector(-wT * (nom.u[end][1] - ℓdes), 0.0, 0.0), (T, 0.0))
    dj = -adj.u[end][3]
    grad = Dict{String, Any}[]
    for ε in (1e-2, 1e-3, 1e-4)
        fd = (objetivo_semilla(ε) - objetivo_semilla(-ε)) / (2ε)
        push!(grad, Dict{String, Any}("eps" => ε, "diferencia_centrada" => fd, "error" => abs(fd - dj)))
    end
    res["gradiente_adjunto"] = Dict{String, Any}("derivada" => dj, "diferencias" => grad)
    @assert grad[end]["error"] < 1e-7

    # Gramianos de controlabilidad y de observabilidad (salida P) sobre una trayectoria nominal.
    function gramianos(inicial, ufun, n)
        pg = con(p; n)
        nominal = integrar_semilla((x, t) -> f(x, ufun(t), pg), inicial, (0.0, T))
        Cy = SMatrix{1, 2}(0.0, 1.0)
        function campo(y, t)
            x = nominal(t)
            A, B = ∂f(x, ufun(t), pg), Bu(x, pg)
            X, W = SMatrix{2, 2}(y[1], y[2], y[3], y[4]), SMatrix{2, 2}(y[5], y[6], y[7], y[8])
            vcat(vec(A * X), vec(A * W + W * A' + B * B'), vec(X' * Cy' * Cy * X))
        end
        y = integrar_semilla(campo, vcat(SVector(1.0, 0.0, 0.0, 1.0), zeros(SVector{8})), (0.0, T)).u[end]
        Dict{String, Any}("controlabilidad" => eigvals(Symmetric([y[5] y[7]; y[6] y[8]])),
                          "observabilidad" => eigvals(Symmetric([y[9] y[11]; y[10] y[12]])))
    end
    res["gramianos_trayectoria_presurizada"] = gramianos(x0, control, 2.0)
    res["gramianos_reposo_n2"] = gramianos(SVector(1.0, 0.0), _ -> 0.0, 2.0)
    res["gramianos_reposo_n1"] = gramianos(SVector(1.0, 0.0), _ -> 0.0, 1.0)
    @assert minimum(res["gramianos_trayectoria_presurizada"]["controlabilidad"]) > 0
    @assert minimum(res["gramianos_trayectoria_presurizada"]["observabilidad"]) > 0
    @assert res["gramianos_reposo_n2"]["controlabilidad"][1] == 0
    @assert res["gramianos_reposo_n1"]["observabilidad"][1] == 0

    # Contraejemplo a la invariancia de P ≥ 0 con controles reales sin restricciones.
    neg = integrar_semilla((x, _) -> f(x, -1.0, con(p; n = 1.0)), SVector(1.0, 0.0), (0.0, 0.1))
    res["control_negativo_n1"] = Dict{String, Any}("q" => -1.0, "tiempo" => 0.1, "estado_final" => collect(neg.u[end]))
    @assert neg.u[end][2] < 0 && neg.u[end][1] < 1

    res["alcance"] = "Comprobaciones de identidades y casos sintéticos. Sin calibración, optimización completa, HJB " *
                     "ni certificación global."
    res
end
