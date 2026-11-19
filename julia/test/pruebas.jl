@testset "pruebas" begin

# ### modelo
@testset "jacobiano analítico contra paso complejo" begin
    rng = Xoshiro(4703)
    for n in (1.0, 2.0, 3.5), _ in 1:10
        p = con(PAR; n)
        x = SVector{2}(0.3 .+ 1.2 .* rand(rng, 2))
        u = 0.1 + 0.9 * rand(rng)
        num = hcat([imag.(f(x .+ 1e-25im .* SVector(i == 1, i == 2), u, p)) ./ 1e-25 for i in 1:2]...)
        @test maximum(abs.(∂f(x, u, p) .- num)) < 1e-11
    end
end

@testset "raiz: Illinois en un intervalo con cambio de signo" begin
    @test raiz(x -> x^3 - 2, 0.0, 2.0) ≈ cbrt(2) atol = 1e-15
    @test raiz(x -> cos(x) - x, 0.0, 1.0) ≈ 0.7390851332151607 atol = 1e-15
    @test_throws ArgumentError raiz(x -> x^2 + 1, -1.0, 1.0)
end
@testset "newton2: sistema 2×2 con raíz conocida" begin
    r = newton2(z -> SVector(z[1]^2 + z[2]^2 - 1, z[1] - z[2]), [0.8, 0.5])
    @test r.retcode == "Success"
    @test r.u ≈ [sqrt(0.5), sqrt(0.5)] atol = 1e-13
end

# ### test_dinamica.py
@testset "test_dinamica.py::test_formula_implicita_sin_control_ni_filtracion_n1" begin
    par = con(PAR; n = 1.0, η = 0.0, ℓ0 = 1.0, P0 = 0.3)
    r = simular(par, 0.0, 5.0, 201)
    S0 = (par.c + par.ℓ0^2) * par.P0
    ell = r.ℓ
    lado_izq = 2 * par.c .* (sqrt.(ell) .- 1) .+ 0.4 .* (ell .^ 2.5 .- 1)
    @test ell[end] > 1.5
    @test maximum(abs.(lado_izq .- S0 .* r.t)) < 1e-8
end

@testset "test_dinamica.py::test_balance_de_volumen_sin_filtracion[$n]" for n in [1.0, 2.0, 3.0]
    par = con(PAR; n = n, η = 0.0)
    U = 0.7
    r = simular(par, U, 3.0, 101)
    S0 = (par.c + par.ℓ0^2) * par.P0
    @test maximum(abs.(r.S .- (S0 .+ U .* r.t))) < 1e-9
end

@testset "test_dinamica.py::test_balance_de_volumen_con_filtracion" begin
    par = con(PAR; n = 2.0, η = 0.5)
    U = 0.5
    r = simular(par, U, 3.0, 101)
    S0 = (par.c + par.ℓ0^2) * par.P0
    @test r.I[end] > 0.1
    @test maximum(abs.(r.S .- S0 .- (U .* r.t .- par.η .* r.I))) < 1e-9
end

@testset "test_dinamica.py::test_barrido_monotonia_positividad_y_cota[$U-$eta-$n]" for U in [0.0, 0.5, 1.0], eta in [0.0, 0.2, 1.0], n in [1.0, 2.0, 3.0]
    par = con(PAR; n = n, η = eta)
    r = simular(par, U, 10.0, 401)
    S0 = (par.c + par.ℓ0^2) * par.P0
    @test length(r.t) == 401 && isapprox(r.t[end], 10.0; rtol = 1e-6, atol = 1e-12)
    @test minimum(diff(r.ℓ)) >= -1e-12
    @test minimum(r.P) >= -1e-12
    @test maximum(r.P .- (S0 .+ U .* r.t) ./ par.c) <= 1e-12
end

@testset "test_dinamica.py::test_equilibrio_sin_presion_ni_control[$eta-$n]" for eta in [0.0, 0.2], n in [1.0, 2.0, 3.0]
    par = con(PAR; n = n, η = eta, P0 = 0.0)
    r = simular(par, 0.0, 5.0, 51)
    @test maximum(abs.(r.ℓ .- par.ℓ0)) < 1e-14
    @test maximum(abs.(r.P)) < 1e-14
end

@testset "test_dinamica.py::test_pendientes_tardias_loglog[$n]" for n in [1.0, 2.0, 3.0]
    par = con(PAR; n = n, η = 0.0)
    T = 1e4
    r = simular(par, 1.0, T, 11)
    pend_ell = log(r.ℓ[end] / r.ℓ[2]) / log(10.0)
    pend_P = log(r.P[end] / r.P[2]) / log(10.0)
    @test abs(pend_ell - 2 * (n + 1) / (3 * n + 2)) <= 0.05
    @test abs(pend_P + (n + 2) / (3 * n + 2)) <= 0.05
end

# ### §3.1 · el largo límite sin inyección, contra la derivación cerrada del informe
# §3.1 deriva que con q ≡ 0 la fractura se detiene donde S se anula, y que ese largo límite es la única
# raíz de n·η·∫ s^(1-n/2)(c+s²)^(n-1) ds = S₀^n, con la cota explícita [ℓ₀^(3n/2) + (3/2η)S₀^n]^(2/(3n))
# y con IGUALDAD en n = 1. Acá se cotejan las tres afirmaciones contra el integrador del proyecto.
# Sin esto, el argumento del informe solo estaba verificado por una corrida suelta de scipy, fuera del
# arnés. Tolerancias MEDIDAS, no elegidas: cinco de las seis combinaciones coinciden a 1e-16 y la lenta
# (n = 1, η = 0,05, donde ℓ_∞ = 7,12) a 3,4e-7, porque a τ = 2000 todavía le queda P = 6e-9.
# `ell_limite` y `cota_ell_limite` viven desde el 15-nov en el paquete (src/objetivo1.jl), porque la
# figura F0 marca esa raíz y esa cota en el plano de fase: con una copia acá, la prueba y la figura
# podían separarse sin que nada avisara.

@testset "§3.1: el largo límite es la raíz de la integral, y la cota lo acota [$n-$eta]" for n in [1.0, 2.0, 3.0], eta in [0.05, 0.2]
    par = con(PAR; n = n, η = eta)
    ell_inf, S0 = ell_limite(par)
    cota = (par.ℓ0^(3n / 2) + 3 / (2eta) * S0^n)^(2 / (3n))
    r = simular(par, 0.0, 2000.0, 11)
    @test abs(r.ℓ[end] - ell_inf) / ell_inf < 1e-5     # el integrador se estanca en la raíz
    @test r.P[end] < 1e-6                              # y se estanca porque S se anuló, no por el horizonte
    @test ell_inf <= cota * (1 + 1e-12)                # la cota explícita
    if n == 1.0
        @test abs(ell_inf - cota) / cota < 1e-10       # con igualdad en n = 1, como dice la derivación
    end
end

@testset "§3.1: con los valores nominales no se llega al objetivo sin inyectar" begin
    # Es lo que hace no trivial el problema de control del §6: la afirmación vale para el caso base,
    # NO en general — con n = 1 el largo límite se va muy por encima de ℓd (3,12 y 7,12 arriba).
    ell_inf, _ = ell_limite(PAR)
    @test ell_inf < PAR.ℓd
    @test ell_limite(con(PAR; n = 1.0, η = 0.2))[1] > PAR.ℓd
end

# ### §3.1 · F0: la figura del plano de fase, contra lo que §3.1 deriva
# La figura del objetivo 1 no sale de un json (el del objetivo 1 es del avance, congelado): sus cifras
# las escribe el generador. Esto es lo que las verifica. Tolerancias MEDIDAS el 15-nov sobre las cuatro
# trayectorias libres: la peor distancia relativa al equilibrio es 6,2e-12 (P₀ = 1) y la peor presión
# sobrante a τ = 200 es 8,7e-7, así que 1e-9 y 1e-5 dejan 160x y 11x de margen.
@testset "§3.1 · F0: la figura del plano de fase dibuja lo que §3.1 deriva" begin
    d = datos_figura_obj1()
    @test length(d.libres) == 4 && length(d.forzadas) == 4
    @testset "familia libre P0 = $(c.P0)" for c in d.libres
        @test abs(c.ell_fin - c.ell_inf) / c.ell_inf < 1e-9   # el integrador llegó al equilibrio
        @test c.P_fin < 1e-5                                   # y llegó porque S se anuló
        @test c.ell_inf <= c.cota * (1 + 1e-12)                # la cota (3.5) acota la raíz (3.4)
        # La cota que la figura MARCA es la misma fórmula que el informe deriva, transcrita aparte.
        q = con(PAR; P0 = c.P0)
        @test cota_ell_limite(q, c.S0) == (q.ℓ0^(3q.n / 2) + 3 / (2q.η) * c.S0^q.n)^(2 / (3q.n))
        @test c.ℓ[1] == PAR.ℓ0 && c.P[1] == c.P0
        @test minimum(diff(c.ℓ)) ≥ -1e-12                      # ℓ no decrece: §3.1, positividad de F
    end
    # El caso base es el que hace no trivial el §6: sin inyectar no se llega a ℓd.
    base = d.libres[2]
    @test base.P0 == PAR.P0
    @test base.ell_inf < PAR.ℓd
    # ℓ_∞ crece con el agua inicial, que es la dependencia S₀^(2/3) que §3.1 lee de la cota.
    @test issorted([c.ell_inf for c in d.libres])
    # El umbral que la figura imprime: ℓ(T; U*) = ℓd, y es único porque ℓ(T; ·) crece.
    @test abs(d.ellT_Ustar - PAR.ℓd) < 1e-10
    @test 0 < d.Ustar < PAR.qmax
    @test all(diff([simular(PAR, U, PAR.T, 3).ℓ[end] for U in range(0.0, PAR.qmax, length = 81)]) .> 0)
    @test all(c.ellT < PAR.ℓd for c in d.forzadas if c.U < d.Ustar)
    @test all(c.ellT > PAR.ℓd for c in d.forzadas if c.U > d.Ustar)
    @test any(c.U < d.Ustar for c in d.forzadas) && any(c.U > d.Ustar for c in d.forzadas)
    # Ninguna trayectoria forzada se detiene: con U > 0 no hay equilibrio (§3.1).
    @test all(c.PT > 0 for c in d.forzadas)
end

# ### test_directo.py
@testset "test_directo.py::test_forward_coincide_con_integracion_de_un_tramo" begin
    q = fill(0.5, 20)
    J, _ = avance(q, PAR)
    s = integrar((z, _) -> push(f(z, 0.5, PAR), L(z, 0.5, PAR)), SVector(PAR.ℓ0, PAR.P0, 0.0), (0.0, PAR.T))
    yT = s.u[end]
    J_ref = Φ(yT[1:2], PAR) + yT[3]
    @test abs(J - J_ref) < 1e-9
end

@testset "test_directo.py::test_gradiente_adjunto_contra_diferencias_finitas" begin
    q = 0.4 .+ 0.2 .* Z4703_A
    J, g = costo_y_gradiente(q, PAR)
    eps = 1e-5
    for k in [1, 6, 12, 16]
        e = zeros(length(q))
        e[k] = eps
        fd = (avance(q + e, PAR)[1] - avance(q - e, PAR)[1]) / (2 * eps)
        @test abs(fd - g[k]) < 1e-7 * max(1.0, abs(g[k]))
    end
    v = Z4703_B
    fd = (avance(q + eps * v, PAR)[1] - avance(q - eps * v, PAR)[1]) / (2 * eps)
    @test abs(fd - dot(g, v)) < 1e-7 * max(1.0, abs(dot(g, v)))
end

@testset "test_directo.py::test_gradiente_detecta_signo_cambiado" begin
    q = fill(0.5, 8)
    J, g = costo_y_gradiente(q, PAR)
    g_malo = PAR.wq .* q .* (PAR.T / 8) .* 2 .- g
    e = zeros(8)
    e[4] = 1e-5
    fd = (avance(q + e, PAR)[1] - avance(q - e, PAR)[1]) / 2e-5
    @test abs(fd - g[4]) < 1e-7
    @test abs(fd - g_malo[4]) > 1e-3
end

@testset "test_directo.py::test_optimo_baja_el_costo_y_anula_el_gradiente" begin
    res = optimizar(PAR, 20)
    d = diagnosticos(res.x, PAR)
    @test res.exito
    @test d["norma_grad_inf"] < 1e-8
    @test d["J"] < avance(zeros(20), PAR)[1]
    @test d["J"] < avance(fill(0.5, 20), PAR)[1]
end

# ### test_cotas_pref.py
@testset "test_cotas_pref.py::test_params_trae_qmax_y_pref" begin
    @test PAR.qmax == 2.0
    @test PAR.Pref == 0.0
end

@testset "test_cotas_pref.py::test_costo_corriente_usa_pref" begin
    par = con(PAR; Pref = 0.5)
    x, u = [1.3, 0.2], 0.7
    @test isapprox(L(x, u, par), 0.5 * (par.wp * 0.3^2 + par.wq * 0.49); atol = 1e-15, rtol = 0)
end

@testset "test_cotas_pref.py::test_dcosto_dx_es_la_derivada_del_costo" begin
    par = con(PAR; Pref = 0.5)
    x, u, h = [1.3, 0.2], 0.7, 1e-6
    fd = [(L(x + h * e, u, par) - L(x - h * e, u, par)) / (2 * h) for e in eachcol(Matrix(1.0I, 2, 2))]
    @test allclose(∂L(x, par), fd; atol = 1e-9)
    @test allclose(∂L(x, par), [0.0, par.wp * (0.2 - 0.5)]; atol = 1e-15)
end

@testset "test_cotas_pref.py::test_gradiente_con_pref_contra_diferencias_finitas" begin
    par = con(PAR; Pref = 0.5)
    q = 0.4 .+ 0.2 .* Z4703_A
    _, g = costo_y_gradiente(q, par)
    eps = 1e-5
    for k in [1, 6, 12, 16]
        e = zeros(length(q))
        e[k] = eps
        fd = (avance(q + e, par)[1] - avance(q - e, par)[1]) / (2 * eps)
        @test abs(fd - g[k]) < 1e-7 * max(1.0, abs(g[k]))
    end
end

@testset "test_cotas_pref.py::test_optimize_respeta_cotas_y_las_activa_con_qmax_08" begin
    par = con(PAR; qmax = 0.8)
    res = optimizar(par, 25)
    @test all(res.x .>= 0.0) && all(res.x .<= 0.8)
    d = diagnosticos(res.x, par)
    @test d["n_cota_sup"] >= 1
    @test d["pmp_brecha_max"] < 5e-3
end

@testset "test_cotas_pref.py::test_diagnostics_cuenta_cotas_en_caso_base" begin
    res = optimizar(PAR, 25)
    d = diagnosticos(res.x, PAR)
    @test d["n_cota_sup"] == 0
    @test d["n_cota_inf"] == 0
end

@testset "test_cotas_pref.py::test_regresion_J200_con_qmax_2" begin
    res, _ = cascada(PAR, (25, 50, 100, 200))
    @test abs(res.fun - 0.37186254275) <= 1e-9
end

# ### test_disparo.py
@testset "test_disparo.py::test_u_estrella_es_el_clip" begin
    x = [1.2, 0.3]
    Dn = PAR.c + 1.44
    for (zP, esperado) in [(0.1, 0.1 / (PAR.wq * Dn)), (-1.0, 0.0), (50.0, PAR.qmax)]
        @test abs(ustar(x, [0.0, zP], PAR) - esperado) < 1e-15
    end
end

@testset "test_disparo.py::test_disparo_penalizado_caso_base" begin
    r = resolver(PAR, SEMILLA, :penalizado)
    @test r["exito"]
    @test r["norma_residuo"] < 1e-10
    @test abs(r["J"] - 0.37186200) <= 1e-8
    for (N, J) in J_DIRECTO
        @test r["J"] < J
    end
    @test maximum(abs.(r["z0"] .- SEMILLA)) <= 5e-3
    @test r["H_oscilacion"] <= 1e-8
end

@testset "test_disparo.py::test_qmax_08_directo_contra_disparo" begin
    par = con(PAR; qmax = 0.8)
    r = resolver(par, SEMILLA, :penalizado)
    @test r["exito"] && r["norma_residuo"] < 1e-10
    res, q = cascada(par, (25, 50, 100, 200))
    @test abs(res.fun - r["J"]) <= 1e-5
    h = par.T / 200
    en_cota = findall(q .>= par.qmax - 1e-12)
    @test length(en_cota) > 0 && all(diff(en_cota) .== 1)
    t_directo = en_cota[end] * h
    @test !isempty(r["conmutaciones"])
    @test abs(r["conmutaciones"][1] - t_directo) <= h
end

@testset "test_disparo.py::test_variante_terminal_y_continuacion" begin
    rt = resolver(PAR, SEMILLA, :terminal)
    @test rt["exito"]
    @test abs(rt["ell_T"] - PAR.ℓd) < 1e-10
    z0, errores = SEMILLA, Float64[]
    for wT in (10.0, 1e2, 1e3, 1e4)
        rp = resolver(con(PAR; wT = wT), z0, :penalizado)
        @test rp["exito"] && rp["norma_residuo"] < 1e-10
        z0 = rp["z0"]
        push!(errores, abs(-wT * (rp["ell_T"] - PAR.ℓd) - rt["zeta_ell_T"]))
    end
    razones = errores[1:end-1] ./ errores[2:end]
    @test all((razones .>= 5) .& (razones .<= 20))
end

@testset "test_disparo.py::test_extremal_da_solucion_densa" begin
    e = extremal(SEMILLA, PAR)
    y = reduce(hcat, [e(τ) for τ in [0.0, PAR.T / 2, PAR.T]])
    @test size(y) == (5, 3)
    @test allclose(y[1:2, 1], [PAR.ℓ0, PAR.P0]) && allclose(y[3:4, 1], SEMILLA)
    @test y[5, 1] == 0.0
end

# ### test_comparaciones.py
@testset "test_comparaciones.py::test_factibilidad_constante_por_brentq" begin
    U, ell_T = constante_factible(PAR)
    @test 0.0 <= U <= PAR.qmax
    @test abs(ell_T - PAR.ℓd) < 1e-10
end

@testset "test_comparaciones.py::test_orden_de_costos" begin
    J_opt = resolver(PAR, SEMILLA, :penalizado)["J"]
    rc = mejor_constante(PAR)
    rr = mejor_rampa(PAR)
    rv = constante_igual_volumen(PAR, SEMILLA)
    @test J_opt < rr["J"] <= rc["J"] + 1e-12
    @test J_opt < rv["J"]
    for r in (rc, rr, rv)
        @test 0.0 <= r["U_min"] && r["U_max"] <= PAR.qmax
    end
end

@testset "test_comparaciones.py::test_desglose_suma_J" begin
    q = collect(range(1.2, 0.1, length = 200))
    d = desglose(q, PAR)
    @test issubset(Set(["terminal", "presion", "control", "volumen", "ell_T", "P_T", "J"]), Set(keys(d)))
    @test abs(d["terminal"] + d["presion"] + d["control"] - d["J"]) < 1e-10
    @test abs(d["control"] - 0.5 * PAR.wq * (sum(q .^ 2) / length(q)) * PAR.T) < 1e-13
    @test abs(d["volumen"] - (sum(q) / length(q)) * PAR.T) < 1e-13
end

@testset "test_comparaciones.py::test_desglose_presion_contra_una_sola_integracion" begin
    par = con(PAR; Pref = 0.25)
    U = 0.6
    s = integrar((z, _) -> push(f(z, U, par), 0.5 * par.wp * (z[2] - par.Pref)^2), SVector(par.ℓ0, par.P0, 0.0),
                 (0.0, par.T))
    yT = s.u[end]
    d = desglose(fill(U, 50), par)
    @test abs(d["presion"] - yT[3]) < 1e-9
    @test abs(d["P_T"] - yT[2]) < 1e-9
    @test abs(d["ell_T"] - yT[1]) < 1e-9
end

@testset "test_comparaciones.py::test_constante_de_igual_volumen_gasta_el_agua_del_optimo" begin
    rv = constante_igual_volumen(PAR, SEMILLA)
    @test abs(rv["U"] * PAR.T - rv["volumen"]) < 1e-14
    _, q = cascada(PAR, (25, 50, 100))
    V_directo = sum(q) * PAR.T / length(q)
    @test abs(rv["volumen"] - V_directo) < 1e-3 * V_directo
end

# ### test_hjb.py
@testset "test_hjb.py::test_interpolacion_bilineal_exacta_en_el_interior" begin
    ejes = _ejes(-1.0, 2.0, 13)
    G0 = [x for x in ejes[1], y in ejes[2]]
    G1 = [y for x in ejes[1], y in ejes[2]]
    V = 0.3 .+ 1.7 .* G0 .- 0.4 .* G1 .+ 2.2 .* G0 .* G1
    rng = Xoshiro(0)
    X0 = -1 .+ 3 .* rand(rng, 500)
    X1 = -1 .+ 3 .* rand(rng, 500)
    exacto = 0.3 .+ 1.7 .* X0 .- 0.4 .* X1 .+ 2.2 .* X0 .* X1
    got = [interpolar(V, ejes, X0[i], X1[i], (:extrapolar, :extrapolar)) for i in 1:500]
    @test maximum(abs.(got .- exacto)) < 1e-12
end

@testset "test_hjb.py::test_extrapolar_es_lineal_y_recortar_congela" begin
    ejes = _ejes(0.0, 1.0, 11)
    G0 = [x for x in ejes[1], y in ejes[2]]
    G1 = [y for x in ejes[1], y in ejes[2]]
    V = 2.0 .* G0 .+ 3.0 .* G1
    X0 = [-0.5, 1.4, 1.9]
    X1 = [-0.3, 1.2, 0.5]
    got = [interpolar(V, ejes, X0[i], X1[i], (:extrapolar, :recortar)) for i in 1:3]
    esperado = 2.0 .* X0 .+ 3.0 .* clamp.(X1, 0.0, 1.0)
    @test maximum(abs.(got .- esperado)) < 1e-12
end

@testset "test_hjb.py::test_sin_dinamica_con_costo_constante" begin
    ejes = _ejes(0.0, 1.0, 11)
    f0 = (X, U) -> (0.0, 0.0)
    costo = (X, U) -> 1.0
    g = X -> X[1]
    V, pol = semi_lagrangiano(f0, costo, g, ejes, [0.0, 1.0], 2.0, 8, (:extrapolar, :recortar))
    G0 = [x for x in ejes[1], y in ejes[2]]
    for k in 0:8
        @test maximum(abs.(V[k+1] .- (G0 .+ 2.0 .- k * 0.25))) < 1e-12
    end
    @test (length(pol), size(pol[1])...) == (8, 11, 11)
end

@testset "test_hjb.py::test_lq_doble_integrador_contra_riccati_baja_al_refinar" begin
    e_grueso = _error_lq(41, 10, 41)
    e_fino = _error_lq(81, 20, 81)
    @test isfinite(e_grueso) && isfinite(e_fino)
    @test e_fino < e_grueso
    @test e_fino < 0.05
end

@testset "test_hjb.py::test_feedback_reproduce_la_politica_en_los_nodos" begin
    ejes = _ejes(-3.0, 3.0, 21)
    sol = solucion(_lq_f, _lq_costo, _lq_g, ejes, range(-4.0, 4.0, length = 21), T_LQ, 5, (:extrapolar, :extrapolar))
    for k in (0, 2, 4)
        for (i, j) in ((11, 11), (8, 14), (4, 17))
            u = feedback(sol, k * sol.dt, (ejes[1][i], ejes[2][j]))
            @test u == sol.politica[k+1][i, j]
        end
    end
end

@testset "test_hjb.py::test_feedback_entre_nodos_propaga_el_tiempo_restante" begin
    ejes = _ejes(-3.0, 3.0, 21)
    sol = solucion(_lq_f, _lq_costo, _lq_g, ejes, range(-4.0, 4.0, length = 21), T_LQ, 5, (:extrapolar, :extrapolar))
    # El control que elige la tabla de Bellman al propagar un horizonte h desde x con V_{k+1}; feedback usa h = t_{k+1} - t.
    function elegir(k, x, h)
        mejor, umejor = Inf, sol.controles[1]
        for u in sol.controles
            e = tabla(sol.f, sol.costo, x, u, h, sol.ejes, sol.bordes)
            q = e[5] + evaluar(sol.V[k+2], e)
            if q < mejor
                mejor, umejor = q, u
            end
        end
        umejor
    end
    casos = [(k, (k + a) * sol.dt, (ejes[1][i], ejes[2][j])) for k in 0:sol.nt-1, a in (0.25, 0.5, 0.75), i in 3:3:19, j in 3:3:19]
    difieren = [(k, t, x) for (k, t, x) in casos if elegir(k, x, (k + 1) * sol.dt - t) != elegir(k, x, sol.dt)]
    @test !isempty(difieren)
    @test all(feedback(sol, t, x) == elegir(k, x, (k + 1) * sol.dt - t) for (k, t, x) in difieren)
    @test all(feedback(sol, t, x) != elegir(k, x, sol.dt) for (k, t, x) in difieren)
end

@testset "test_hjb.py::test_feedback_solo_acepta_t_en_0_T" begin
    ejes = _ejes(-3.0, 3.0, 21)
    sol = solucion(_lq_f, _lq_costo, _lq_g, ejes, range(-4.0, 4.0, length = 21), T_LQ, 5, (:extrapolar, :extrapolar))
    x = (ejes[1][11], ejes[2][11])
    @test_throws ArgumentError feedback(sol, sol.T, x)
    @test_throws ArgumentError feedback(sol, sol.T + 0.1, x)
    @test_throws ArgumentError feedback(sol, -0.1, x)
    @test feedback(sol, 0.0, x) isa Float64
    @test feedback(sol, prevfloat(sol.T), x) isa Float64
end

@testset "test_hjb.py::test_lazo_con_control_constante_contra_solucion_exacta" begin
    u, x0, T = 0.7, [1.0, -0.5], 1.0
    Ainv = inv(A_OSC)
    x_exacto(t) = (E = exp(A_OSC * t); E * x0 + Ainv * (E - I) * B_OSC * u)
    J_exacto = quadgk(t -> _osc_costo(x_exacto(t), u), 0.0, T; atol = 1e-14, rtol = 1e-13)[1] + _osc_g(x_exacto(T))
    r = simular_lazo(_osc_f, _osc_costo, _osc_g, (t, x) -> u, x0, T; n_sub = 40, n_rk = 5)
    @test maximum(abs.(r.X[end] .- x_exacto(T))) < 1e-10
    @test abs(r.J - J_exacto) < 1e-10
end

@testset "test_hjb.py::test_lazo_con_la_ley_lq_exacta_contra_riccati" begin
    x0, T = [1.0, -0.5], 1.0
    ric = _riccati_osc(T)
    J_star = 0.5 * dot(x0, reshape(ric(0.0), 2, 2) * x0)
    ley(t, x) = float(-(B_OSC' * reshape(ric(t), 2, 2) * x))
    e = [simular_lazo(_osc_f, _osc_costo, _osc_g, ley, x0, T; n_sub = n).J - J_star for n in (200, 400)]
    @test e[1] > 0 && e[2] > 0
    @test e[2] < 1e-6
    @test 3.0 < e[1] / e[2] < 5.0
end

@testset "test_hjb.py::test_hjb_fractura_dominio_y_valor_en_x0" begin
    par = Par()
    rd = resolver(par, [0.8161, 0.4346], :penalizado)
    Y = reduce(hcat, [rd["sol"](τ) for τ in range(0.0, par.T, length = 2001)])
    sol = hjb_fractura(par, 81, 81, 20, 41)
    (l0, l1), (p0, p1) = [(e[1], e[end]) for e in sol.ejes]
    @test l0 <= par.ℓ0 && l1 >= 1.2 * maximum(Y[1, :])
    @test p0 <= minimum(Y[2, :]) && p1 >= 1.5 * maximum(Y[2, :])
    i = argmin(abs.(sol.ejes[1] .- par.ℓ0))
    j = argmin(abs.(sol.ejes[2] .- par.P0))
    V0 = sol.V[1][i, j]
    @test abs(V0 - rd["J"]) < 0.01 * rd["J"]
    lc = lazo_cerrado(par, sol)
    @test 0.0 <= lc.J - rd["J"] < 5e-3
    @test lc.J - rd["J"] < 5e-4   # H01: propagando un dt completo daba 2,17e-3; con el tiempo restante, 1,59e-4
end

# ### objetivo2.jl · F1 tareas 1.1 a 1.3 (sistema linealizado)
# Criterios fijados ANTES de correr, en el plan 2026-11-14 (F1). La nominal de referencia para el
# control cruzado es la MISMA de la semilla: p = Par(), T = 2, x0 = (1, 0.3), u(t) = 0.5 + 0.1 sin t,
# más las dos de reposo (x = (1,0), u ≡ 0) con n = 2 y n = 1, que son las tres que la semilla calcula.
U_SEM(t) = 0.5 + 0.1 * sin(t)
X0_SEM = SVector(1.0, 0.3)
_nom_sem() = nominal_control(con(PAR; n = 2.0), U_SEM, X0_SEM; T = 2.0)
_nom_reposo(n) = nominal_control(con(PAR; n), _ -> 0.0, SVector(1.0, 0.0); T = 2.0)
# Memoizada: cada `semilla_lectura()` recorre las 8 identidades en BigFloat, y seis testsets la piden.
_SEM_CACHE = Ref{Any}(nothing)
_sem() = (_SEM_CACHE[] === nothing && (_SEM_CACHE[] = semilla_lectura()); _SEM_CACHE[])

@testset "1.1 A(t) y B(t) sobre la nominal contra diferencias finitas" begin
    nom = _nom_sem()
    for t in range(0.0, 2.0, length = 20)
        A, B = matrices_lineales(nom, t)
        x, u = nom.x(t), nom.u(t)
        h = 1e-6
        Anum = hcat(((f(x + h * e, u, nom.p) - f(x - h * e, u, nom.p)) / 2h
                     for e in (SVector(1.0, 0.0), SVector(0.0, 1.0)))...)
        Bnum = (f(x, u + h, nom.p) - f(x, u - h, nom.p)) / 2h
        @test maximum(abs.(A - Anum)) < 1e-7 * max(1.0, maximum(abs, A))
        @test maximum(abs.(B - Bnum)) < 1e-7 * max(1.0, maximum(abs, B))
    end
end

@testset "1.1 Φ: propiedad de grupo y Liouville" begin
    nom = _nom_sem()
    t0, t1, t2 = 0.0, 0.7, 1.9
    @test maximum(abs.(fundamental(nom, t1, t2) * fundamental(nom, t0, t1) -
                       fundamental(nom, t0, t2))) < 1e-8
    # det Φ(t,0) = exp(∫₀ᵗ tr A): control independiente del integrador de matrices.
    for t in (0.5, 1.0, 2.0)
        tr∫ = quadgk(s -> tr(matrices_lineales(nom, s)[1]), 0.0, t; rtol = 1e-12, atol = 1e-14)[1]
        @test abs(det(fundamental(nom, 0.0, t)) - exp(tr∫)) < 1e-8 * exp(tr∫)
    end
    # Ψ(s) = Φ(T,s) integrado hacia atrás: Ψ(s)Φ(s,0) = Φ(T,0), sin invertir nada.
    Ψ, Φ = propagador_final(nom, 2.0), fundamental_densa(nom, 0.0, 2.0)
    ΦT = Φ(2.0)
    for s in (0.0, 0.3, 1.1, 2.0)
        @test maximum(abs.(Ψ(s) * Φ(s) - ΦT)) < 1e-8 * max(1.0, maximum(abs, ΦT))
    end
end

@testset "1.2 Gramiano de controlabilidad contra la semilla ($caso)" for caso in
        ("trayectoria_presurizada", "reposo_n2", "reposo_n1")
    nom = caso == "trayectoria_presurizada" ? _nom_sem() :
          _nom_reposo(caso == "reposo_n1" ? 1.0 : 2.0)
    W = gramiano_controlabilidad(nom; T = 2.0)
    @test maximum(abs.(W - W')) < 1e-12 * max(1.0, maximum(abs, W))
    ref = _sem()["gramianos_$caso"]["controlabilidad"]          # valores propios, BigFloat → Float64
    mios = eigvals(Symmetric(Matrix(W)))
    # La semilla integra la EDO de Lyapunov W' = AW + WA' + BB'; acá es cuadratura sobre Φ(T,s).
    # Dos caminos distintos: el criterio es 1e-8 relativo al mayor valor propio.
    @test maximum(abs.(mios - ref)) < 1e-8 * max(1.0, maximum(abs, ref))
    @test minimum(mios) > -1e-12 * max(1.0, maximum(mios))      # semidefinida positiva
end

@testset "1.2 energía mínima: alcanza el blanco y su energía es zᵀW⁻¹z" begin
    nom = _nom_sem()
    W = gramiano_controlabilidad(nom; T = 2.0)
    zT = SVector(0.013, -0.021)
    em = energia_minima(nom, zT; T = 2.0)
    # Integrar el sistema LINEAL con ese control: z(T) tiene que dar zT.
    z = integrar((y, t) -> let (A, B) = matrices_lineales(nom, t); A * y + B * em.u(t) end,
                 SVector(0.0, 0.0), (0.0, 2.0)).u[end]
    @test maximum(abs.(z - zT)) < 1e-8 * max(1.0, maximum(abs, zT))
    @test abs(em.energia - dot(zT, Symmetric(Matrix(W)) \ Vector(zT))) < 1e-9 * abs(em.energia)
    # Y su energía es la mínima: ∫u² = zᵀW⁻¹z, medida por cuadratura sobre el control devuelto.
    @test abs(quadgk(s -> em.u(s)^2, 0.0, 2.0; rtol = 1e-12)[1] - em.energia) < 1e-8 * abs(em.energia)
end

@testset "1.2 Kalman congelado: la degeneración vive en P = 0" begin
    # Criterio fijado antes: en reposo (P ≡ 0) con n = 2, ∂F/∂ℓ = ∂F/∂P = 0 y la primera fila de A
    # se anula, así que [B AB] pierde rango. Fuera de P = 0 el rango es 2.
    k0 = kalman_congelado(_nom_reposo(2.0), 1.0)
    @test k0.rango == 1 && k0.sv[2] / k0.sv[1] < 1e-10
    kp = kalman_congelado(_nom_sem(), 1.0)
    @test kp.rango == 2 && kp.sv[2] / kp.sv[1] > 1e-10
end

@testset "1.3 Gramiano de observabilidad contra la semilla ($caso)" for caso in
        ("trayectoria_presurizada", "reposo_n2", "reposo_n1")
    nom = caso == "trayectoria_presurizada" ? _nom_sem() :
          _nom_reposo(caso == "reposo_n1" ? 1.0 : 2.0)
    M = gramiano_observabilidad(nom; T = 2.0)
    @test maximum(abs.(M - M')) < 1e-12 * max(1.0, maximum(abs, M))
    ref = _sem()["gramianos_$caso"]["observabilidad"]
    @test maximum(abs.(eigvals(Symmetric(Matrix(M))) - ref)) < 1e-8 * max(1.0, maximum(abs, ref))
end

@testset "1.3 la presión observa ℓ: norma L² = forma cuadrática del Gramiano" begin
    nom = _nom_sem()
    M = gramiano_observabilidad(nom; T = 2.0)
    δ = SVector(1.0, 0.0)                     # perturbación solo en ℓ
    # (a) En el sistema lineal es una identidad: ∫|CΦδ|² ds = δᵀMδ. Control de implementación.
    Φ = fundamental_densa(nom, 0.0, 2.0)
    l2 = quadgk(s -> (Φ(s) * δ)[2]^2, 0.0, 2.0; rtol = 1e-12, atol = 1e-14)[1]
    @test abs(l2 - dot(δ, M * δ)) < 1e-8 * abs(l2)
    # (b) En el sistema NO lineal, con ε = 1e-4, la razón tiende a 1 con error O(ε). Esto es lo que
    #     hace honesto el enunciado «se observa ℓ a través de P»: criterio fijado antes, < 1e-2.
    ε = 1e-4
    xa = integrar((y, t) -> f(y, nom.u(t), nom.p), X0_SEM, (0.0, 2.0))
    xb = integrar((y, t) -> f(y, nom.u(t), nom.p), X0_SEM + ε * δ, (0.0, 2.0))
    l2nl = quadgk(s -> ((xb(s)[2] - xa(s)[2]) / ε)^2, 0.0, 2.0; rtol = 1e-10)[1]
    @test abs(l2nl / dot(δ, M * δ) - 1) < 1e-2
    @test l2nl > 0                            # la presión SÍ cambia cuando cambia ℓ
end

@testset "1.3 observabilidad congelada: en reposo ℓ no se observa" begin
    o0 = observabilidad_congelada(_nom_reposo(2.0), 1.0)
    @test o0.rango == 1 && o0.sv[2] / o0.sv[1] < 1e-10
    op = observabilidad_congelada(_nom_sem(), 1.0)
    @test op.rango == 2 && op.sv[2] / op.sv[1] > 1e-10
end

# ### objetivo2.jl · F1 tarea 1.4 (LQ de seguimiento por Riccati)
# Honestidad del ciclo: `riccati` y `lazo_cerrado` se escribieron junto con el resto de objetivo2.jl,
# así que ESTAS pruebas no vieron la fase roja de «función no definida» como las de 1.1–1.3. Lo que sí
# cubren es el error que de verdad se cuela acá: el signo de K, que sin la prueba de perturbación nula
# pasa desapercibido porque el lazo igual converge a algo parecido.
#
# Los pesos NO son de adorno: Qf = diag(w_T, 0), Q = diag(0, w_p) y R = w_q son la segunda variación
# del costo original Φ(x(T)) + ∫L, o sea el LQ que sale de este problema (Gauss–Newton: se omite la
# segunda derivada de f contra el adjunto). Para las pruebas estructurales se usa además Q = Qf = I,
# donde P sí es definida positiva y el criterio se puede escribir sin ambigüedad.
_pesos_problema(p) = (Q = SMatrix{2,2}(0.0, 0.0, 0.0, p.wp), Qf = SMatrix{2,2}(p.wT, 0.0, 0.0, 0.0), R = p.wq)

@testset "1.4 Riccati: simétrica, definida positiva y satisface su EDO" begin
    nom = nominal_optima(PAR, optimizar(PAR, 25).x)
    P = riccati(nom; Q = I2, R = 1.0, Qf = I2)
    for t in range(0.0, PAR.T, length = 20)
        Pt = P(t)
        @test maximum(abs.(Pt - Pt')) < 1e-12 * max(1.0, maximum(abs, Pt))
        @test minimum(eigvals(Symmetric(Matrix(Pt)))) > 0          # definida positiva con Q = Qf = I
    end
    # Residual de -Ṗ = PA + A'P - PBR⁻¹B'P + Q por diferencias finitas centradas.
    h = 1e-5
    for t in range(h, PAR.T - h, length = 20)
        A, B = matrices_lineales(nom, t)
        Pt = P(t)
        res = (P(t + h) - P(t - h)) / 2h + (Pt * A + A' * Pt - (Pt * B) * (B' * Pt) + I2)
        @test maximum(abs.(res)) < 1e-8 * max(1.0, maximum(abs, Pt))
    end
end

@testset "1.4 perturbación nula: el lazo cerrado reproduce la nominal (control de signo de K)" begin
    nom = nominal_optima(PAR, optimizar(PAR, 25).x)
    w = _pesos_problema(PAR)
    P = riccati(nom; w.Q, w.R, w.Qf)
    lc = lazo_cerrado(nom, P; w.R)
    @test maximum(abs.(lc.xT - nom.x(PAR.T))) < 1e-9
    for t in range(0.0, PAR.T, length = 21)
        @test maximum(abs.(lc.x(t) - nom.x(t))) < 1e-9
    end
end

@testset "1.4 el experimento: LQ contra lazo abierto bajo perturbación ($cual)" for cual in ("P0", "eta")
    # Criterio fijado ANTES de correr: el LQ reduce la DESVIACIÓN TERMINAL respecto de la nominal,
    # |ℓ(T) − ℓ̄(T)|, en al menos un factor 3 frente al lazo abierto. El plan escribió el criterio sobre
    # |ℓ(T) − ℓ_d|; se cambia a la desviación de la nominal, y la razón se declara acá antes de medir:
    # el óptimo del directo NO lleva ℓ(T) a ℓ_d (hay compromiso con el costo), así que ℓ_d no es lo que
    # este controlador persigue. Un LQ de seguimiento sigue la nominal. Las dos métricas se reportan.
    nom = nominal_optima(PAR, optimizar(PAR, 25).x)
    w = _pesos_problema(PAR)
    P = riccati(nom; w.Q, w.R, w.Qf)
    x0p = cual == "P0" ? SVector(PAR.ℓ0, 1.2 * PAR.P0) : SVector(PAR.ℓ0, PAR.P0)
    pp = cual == "eta" ? con(PAR; η = 1.2 * PAR.η) : PAR
    abierto = lazo_cerrado(nom, P; x0 = x0p, p = pp, w.R, lazo_abierto = true)
    cerrado = lazo_cerrado(nom, P; x0 = x0p, p = pp, w.R)
    ℓT = nom.x(PAR.T)[1]
    da, dc = abs(abierto.xT[1] - ℓT), abs(cerrado.xT[1] - ℓT)
    @test da > 1e-6                       # control: la perturbación SÍ mueve el estado terminal
    @test dc < da / 3
end

# ### F2 · el acuerdo entre los tres métodos
# `acuerdo()` corre el disparo, la cascada del directo hasta N = 400 y los tres niveles de HJB (~15 s).
# Memoizado con un Ref, como `_sem()`: tres testsets lo piden y correrlo tres veces no agrega información.
_ACU_CACHE = Ref{Any}(nothing)
_acu() = (_ACU_CACHE[] === nothing && (_ACU_CACHE[] = acuerdo()); _ACU_CACHE[])

@testset "2.1 acuerdo: los tres métodos dan el mismo J" begin
    a = _acu()[1]
    @test a["criterios"]["disparo_converge"]
    @test a["criterios"]["tres_metodos_coinciden_en_J"]
    # Control de que la tolerancia SÍ discrimina: el mejor control constante (J = 0,396738, del
    # objetivo 3) está 6,7 % arriba del óptimo, muy fuera de la tolerancia de 5e-4 del criterio.
    @test abs(0.3967375482981293 - a["J_disparo"]) / a["J_disparo"] > 5e-4
end

@testset "2.1 acuerdo: el directo converge con orden 2" begin
    a = _acu()[1]
    @test a["criterios"]["brechas_directo_bajan"]
    @test a["criterios"]["orden_directo_cerca_de_2"]
    @test a["criterios"]["brechas_hjb_bajan"]
end

@testset "2.1 acuerdo: coinciden las CURVAS de control, no solo J" begin
    a, d = _acu()
    @test a["criterios"]["sup_dU_directo_baja"]
    @test a["criterios"]["sup_dU_hjb_baja"]
    @test maximum(f["sup_dU"] for f in a["directo"]) < 0.05
    # Control: la misma métrica contra el mejor control constante (U = 0,791) tiene que dar un
    # número grande. Si no, «sup|ΔU| chico» no estaría midiendo acuerdo entre curvas.
    @test maximum(abs.(d.pmp.U .- 0.7913865447989883)) > 0.25
end

# ---- F3 · sensibilidades, calibre e identificabilidad -------------------------------------------
TS_SENS = collect(range(0.1, 2.0, length = 20))
_u_sens(t) = 0.8 + 0.3 * sin(2t)
TS_DIM = collect(range(0.02, 4.0, length = 200))     # la misma malla del .py de MA5306
_q_uno(_) = 1.0

@testset "3.1 sensibilidades: la EDO aumentada contra diferencias finitas centradas ($par)" for par in
        (:n, :c, :η, :P0)
    p = Par()
    S = sens_presion(p, _u_sens, TS_SENS; pars = (par,))[:, 1]
    θ = Float64(getfield(p, par))
    h = 1e-4 * max(1.0, abs(θ))
    Pmas = presion_adim(con(p; (par => θ + h,)...), _u_sens, TS_SENS)
    Pmen = presion_adim(con(p; (par => θ - h,)...), _u_sens, TS_SENS)
    num = (Pmas .- Pmen) ./ 2h
    esc = max(1e-3, maximum(abs, num))
    @test maximum(abs, S .- num) < 1e-6 * esc
    # Control de que la escala del criterio discrimina: la sensibilidad no es cero.
    @test maximum(abs, num) > 1e-3
end

@testset "3.1 sensibilidades: las cuatro columnas juntas son las cuatro por separado" begin
    p = Par()
    J = sens_presion(p, _u_sens, TS_SENS)
    for (j, par) in enumerate(PARS_SENS)
        @test maximum(abs, J[:, j] - sens_presion(p, _u_sens, TS_SENS; pars = (par,))[:, 1]) <
              1e-9 * max(1.0, maximum(abs, J[:, j]))
    end
    # Relativa = θ_j por la columna: ∂P/∂log θ_j.
    Jr = sens_presion(p, _u_sens, TS_SENS; relativa = true)
    for (j, par) in enumerate(PARS_SENS)
        @test maximum(abs, Jr[:, j] - getfield(p, par) * J[:, j]) < 1e-12 * max(1.0, maximum(abs, Jr[:, j]))
    end
end

@testset "3.1 calibre: el núcleo de dos direcciones es una simetría EXACTA (n = $nn)" for nn in
        (1.5, 2.0, 3.0)
    ref = presion_dim(DIM_BASE, nn, 0.3, _q_uno, TS_DIM)
    esc = maximum(abs, ref)
    esc > 0.1 || error("el caso base dimensional no presuriza: esc = $esc")
    # Amplitud FINITA, no infinitesimal: a = 2 y b = 3 son factores 2 y 3 en los parámetros. Una
    # dirección «casi nula» se separaría acá; una simetría exacta, no.
    for (sa, sb) in ((0.1, 0.0), (log(2.0), 0.0), (0.0, log(3.0)), (log(2.0), log(3.0)))
        θ = mover_log(mover_log(DIM_BASE, calibre_escala(nn), sa), calibre_vk(nn), sb)
        d = presion_dim(θ, nn, 0.3, _q_uno, TS_DIM) .- ref
        @test maximum(abs, d) < 1e-10 * esc
    end
    # Controles negativos: dos direcciones que NO son del calibre mueven la presión 6 a 11 órdenes
    # de magnitud más. Sin esto, «la diferencia es 1e-14» no diría nada.
    for k in (:C, :E)
        dir = merge(map(_ -> 0.0, DIM_BASE), NamedTuple{(k,)}((1.0,)))
        d = presion_dim(mover_log(DIM_BASE, dir, 0.1), nn, 0.3, _q_uno, TS_DIM) .- ref
        @test maximum(abs, d) > 1e-4 * esc
    end
end

@testset "3.1 calibre: sus cuatro invariantes no se mueven, y un quinto monomio sí (n = $nn)" for nn in
        (1.5, 2.0, 3.0)
    for (sa, sb) in ((log(2.0), 0.0), (0.0, log(3.0)), (log(2.0), log(3.0)))
        θ = mover_log(mover_log(DIM_BASE, calibre_escala(nn), sa), calibre_vk(nn), sb)
        for (nombre, w) in invariantes(nn)
            @test abs(monomio(θ, w) - monomio(DIM_BASE, w)) < 1e-12 * abs(monomio(DIM_BASE, w))
        end
        # Control: E solo (que no es invariante) sí se mueve, salvo en la dirección v₀–K, que no lo toca.
        mE = monomio(θ, merge(map(_ -> 0.0, DIM_BASE), (E = 1.0,)))
        sa == 0.0 || @test abs(mE - DIM_BASE.E) > 0.5 * DIM_BASE.E
    end
end

@testset "3.1 calibre: el Jacobiano dimensional lo tiene por núcleo, y su espectro es el del .py" begin
    nn = 2.0
    J = jacobiano_dim(DIM_BASE, nn, 0.3, _q_uno, TS_DIM)
    sv = svdvals(J) ./ svdvals(J)[1]
    # Referencia: la corrida del 14-nov de 2026-10-16-fisher-ma4703.py (q constante, 6 parámetros
    # libres), que imprime 1.0e+00 8.8e-02 4.2e-04 1.7e-05 2.2e-10 2.3e-11. Acá el integrador es
    # Vern9 y allá Radau: coincidir en los cuatro grandes cruza dos implementaciones.
    for (k, ref) in zip(2:4, (8.8e-2, 4.2e-4, 1.7e-5))
        @test abs(sv[k] - ref) < 0.1 * ref
    end
    # Los dos chicos están en el piso de ruido de las diferencias finitas, no en un valor con sentido.
    @test sv[5] < 1e-7 && sv[6] < 1e-7
    # Y el núcleo numérico cae en el subespacio que generan las DOS direcciones derivadas a mano.
    @test residuo_nucleo(J, (calibre_escala(nn), calibre_vk(nn))) < 1e-3
    # Control: un subespacio equivocado (dos direcciones cualesquiera) da residuo de orden 1.
    @test residuo_nucleo(J, (merge(map(_ -> 0.0, DIM_BASE), (E = 1.0,)),
                             merge(map(_ -> 0.0, DIM_BASE), (C = 1.0,)))) > 0.5
end

@testset "3.1 calibre: está CLAVADO; mover un exponente 1 % ya rompe la invariancia (n = $nn)" for nn in
        (2.0, 3.0)
    ref = presion_dim(DIM_BASE, nn, 0.3, _q_uno, TS_DIM)
    esc = maximum(abs, ref)
    # El exponente de L₀ en la dirección de escala es 1 exacto. Si fuera «aproximadamente 1», la
    # invariancia a amplitud finita (a = 2) se caería: esto mide cuánto, y es el control de que el
    # 1e-14 de los testsets anteriores viene del calibre y no de que el modelo sea insensible a L₀.
    errores = Float64[]
    for δ in (0.0, 0.01, 0.05, 0.2)
        dir = merge(calibre_escala(nn), (L0 = 1.0 + δ,))
        push!(errores, maximum(abs, presion_dim(mover_log(DIM_BASE, dir, log(2.0)), nn, 0.3, _q_uno,
                                                TS_DIM) .- ref) / esc)
    end
    @test errores[1] < 1e-10
    @test all(e > 1e-4 for e in errores[2:end])
    # Control de discriminación: el salto entre el exponente exacto y el 1 % es de varios órdenes.
    @test errores[2] > 1e5 * max(errores[1], 1e-16)
    # Y es monótono en δ: no es un número suelto que pasó por casualidad.
    @test issorted(errores)
end

# ---- F3.2 · la matriz de información de Fisher --------------------------------------------------
TS_FISHER = collect(range(0.2, 2.0, length = 40))

@testset "3.2 Fisher: M es simétrica, ⪰ 0, y con un solo instante tiene rango 1" begin
    p = Par()
    M = fisher(p, _u_sens, TS_FISHER)
    @test issymmetric(Matrix(M))
    ev = eigvals(Symmetric(Matrix(M)))
    @test minimum(ev) > -1e-6 * maximum(ev)
    @test rank(Matrix(M); rtol = 1e-10) == length(PARS_SENS)
    # Con un solo instante M = s sᵗ/σ²: rango 1 y traza = |s|²/σ². Si el ensamblado sumara mal
    # (por ejemplo J Jᵗ en vez de Jᵗ J), esto se cae.
    M1 = fisher(p, _u_sens, [1.0])
    @test rank(Matrix(M1); rtol = 1e-10) == 1
    s = vec(sens_presion(p, _u_sens, [1.0]))
    @test abs(tr(Matrix(M1)) - sum(abs2, s) / SIGMA_FISHER^2) < 1e-6 * tr(Matrix(M1))
    # σ fija la escala y nada más: M(σ/2) = 4 M(σ).
    @test maximum(abs, Matrix(fisher(p, _u_sens, TS_FISHER; σ = SIGMA_FISHER / 2)) .-
                       4 .* Matrix(M)) < 1e-6 * maximum(abs, Matrix(M))
    # La inflación es ≥ 1 por construcción (Cauchy–Schwarz) y acá es > 1: los parámetros se apoyan
    # unos en otros. Si diera 1 en todos, leer 1/M_jj sería lo mismo que perfilar.
    @test all(inflacion(M) .>= 1 - 1e-8)
    @test maximum(inflacion(M)) > 1.5
end

@testset "3.2 Fisher: Cramér–Rao predice la dispersión de mínimos cuadrados (200 réplicas)" begin
    p = Par()
    r = experimento_cramer_rao(p, _u_sens, TS_FISHER)
    @test r.iteraciones_max < 20                      # Gauss–Newton convergió en todas las réplicas
    # «Bien identificado» se define ANTES de mirar: cota relativa menor al 10 % del valor.
    bien = findall(r.cota ./ abs.(r.theta) .< 0.10)
    @test !isempty(bien)
    for j in bien
        @test 1 / 1.5 < r.razon[j] < 1.5              # dentro de un factor 1,5 de la cota
        @test abs(r.sesgo_sobre_cota[j]) < 0.5        # el estimador no está sesgado a esta escala
    end
    # Control de que la cota no es un número cualquiera: con σ diez veces más chico, la cota baja
    # diez veces y la RAZÓN empírica/cota no se mueve (es invariante en σ en el régimen lineal).
    r2 = experimento_cramer_rao(p, _u_sens, TS_FISHER; σ = SIGMA_FISHER / 10, replicas = 50)
    @test maximum(abs, r2.cota .* 10 .- r.cota) < 1e-3 * maximum(r.cota)
    for j in bien
        @test 1 / 1.5 < r2.razon[j] < 1.5
    end
    # Control de mutación: comparar contra 1/sqrt(M_jj) —la varianza con los demás FIJOS, que es la
    # trampa— da un número demasiado chico justo en los parámetros que se apoyan en otros.
    ingenua = 1 ./ sqrt.(diag(r.M))
    @test maximum(r.cota ./ ingenua) > 1.5
end

# ---- F3.3 · la frontera control ↔ información ---------------------------------------------------
Q_DIS = [0.6, 0.9, 1.2, 1.0, 0.8, 0.7, 0.5, 0.4]

@testset "3.3 diseño: el escalón del directo y el integrador de sensibilidades dan la misma presión" begin
    p = Par()
    ts = [0.3, 0.7, 1.1, 1.5, 1.9]
    t = malla(p, length(Q_DIS))
    Pa = presion_adim(p, escalon(Q_DIS, p), ts)
    tramos = avance(Q_DIS, p)[2]
    Pd = [tramos[clamp(searchsortedlast(t, τ), 1, length(Q_DIS))](τ)[2] for τ in ts]
    # J usa RK4 por tramos y M usa Vern9 continuo: son dos integradores del mismo sistema, y el
    # objetivo combinado solo tiene sentido si coinciden. Esto lo mide en vez de suponerlo.
    @test maximum(abs, Pa .- Pd) < 1e-7
    # Control: contra otro control la diferencia es cinco órdenes mayor. Sin esto, comparar algo
    # consigo mismo pasaría igual.
    @test maximum(abs, presion_adim(p, escalon(reverse(Q_DIS), p), ts) .- Pd) > 1e-2
end

@testset "3.3 diseño: en γ = 0 el óptimo del diseño ES el óptimo del método directo" begin
    p = Par()
    ref = optimizar(p, NS_DISENO; q0 = fill(0.5, NS_DISENO))
    r = optimizar_diseno(p, 0.0; q0 = fill(0.5, NS_DISENO))
    @test r.exito
    @test abs(r.J - ref.fun) < 1e-6 * max(1.0, abs(ref.fun))
    @test maximum(abs, r.q .- ref.x) < 1e-4
end

@testset "3.3 diseño: la derivada direccional del gradiente ensamblado concuerda con el objetivo" begin
    p = Par()
    γ = 0.5
    rng = Xoshiro(3303)
    d = randn(rng, length(Q_DIS))
    d ./= norm(d)
    fd = (psi(Q_DIS .+ 1e-5 .* d, p, γ) - psi(Q_DIS .- 1e-5 .* d, p, γ)) / 2e-5
    @test abs(dot(grad_psi(Q_DIS, p, γ), d) - fd) < 1e-5 * max(1.0, abs(fd))
    # Control de mutación: con el término de información fuera del gradiente, la misma comparación
    # falla por tres órdenes. Es lo que prueba que γ, el signo y los índices entran bien.
    @test abs(dot(gradiente(Q_DIS, p), d) - fd) > 1e-2 * max(1.0, abs(fd))
    # Y el término pesa: su gradiente no es ruido numérico.
    @test norm(grad_logdet_fisher(Q_DIS, p)) > 1e-2
    # En γ = 0 el objetivo y el gradiente son exactamente los del directo (sin llamar a Fisher).
    @test psi(Q_DIS, p, 0.0) == avance(Q_DIS, p)[1]
    @test grad_psi(Q_DIS, p, 0.0) == gradiente(Q_DIS, p)
end


# ---- F4.1 · el AT2 en 1D: forward --------------------------------------------------------------

@testset "4.1 AT2: la rama homogénea reproduce la solución cerrada" begin
    p = ParAT2(; L = 1.0, ℓ = 1.0, N = 40, defecto = 0.0)
    εs = [0.2, 0.5, 0.8, 1.0, 1.4] .* eps_c(p)
    sol = resolver_at2(εs .* p.L, p)
    @test sol.kkt < 1e-10
    @test !any(any, sol.tope)
    for (i, ε) in enumerate(εs)
        αa = alpha_homogenea(ε, p)
        σa = sigma_homogenea(ε, p)
        @test maximum(abs, sol.alpha[i] .- αa) < 1e-9
        @test abs(sol.sigma[i] - σa) < 1e-9 * σa
    end
    # El pico de la rama homogénea: ε = ε_c/√3 y σ = (3√3/16)√(E G_c/ℓ). Es LA relación conocida del
    # AT2, y acá se mide contra el solver en vez de citarse.
    pk = resolver_at2([eps_pico_homogenea(p) * p.L], p)
    @test abs(pk.sigma[1] - sigma_pico_homogenea(p)) < 1e-9 * sigma_pico_homogenea(p)
    @test sigma_homogenea(0.9 * eps_pico_homogenea(p), p) < sigma_pico_homogenea(p)
    @test sigma_homogenea(1.1 * eps_pico_homogenea(p), p) < sigma_pico_homogenea(p)
    # Control de mutación: con ε_c² = G_c/(2Eℓ) —el factor 2 que se cuela si c_w se normaliza mal— la
    # misma comparación falla por seis órdenes. Sin esto la fórmula se verificaría sola.
    mal = ParAT2(; L = 1.0, ℓ = 2.0, N = 40)
    @test maximum(abs(sigma_homogenea(ε, mal) - sol.sigma[i]) for (i, ε) in enumerate(εs)) > 1e-3
end

@testset "4.1 AT2: la alternada baja la energía y la irreversibilidad congela el daño" begin
    p = ParAT2(; L = 1.0, ℓ = 0.3, N = 60, defecto = 0.0)
    εc = eps_c(p)
    # Sube a 0,9 ε_c, descarga a 0,7 y 0,5, vuelve a 0,85 —todavía bajo el máximo histórico— y recién
    # a 0,95 pasa de largo. Las etapas 4, 5 y 6 son las que miden la irreversibilidad.
    carga = [0.3, 0.6, 0.9, 0.7, 0.5, 0.85, 0.95] .* εc .* p.L
    sol = resolver_at2(carga, p)
    @test sol.kkt < 1e-12
    @test !any(any, sol.tope)
    for k in 2:length(carga)
        @test minimum(sol.alpha[k] .- sol.alpha[k - 1]) > -1e-14
    end
    # En las tres etapas bajo el máximo histórico α queda congelado BIT A BIT, y el conjunto activo es
    # todo el dominio: ahí la derivada del estado respecto de la carga es cero, y eso es lo que el
    # adjunto de F4.2 congela.
    for k in (4, 5, 6)
        @test all(sol.activos[k])
        @test sol.alpha[k] == sol.alpha[3]
    end
    # Y al pasar de largo se descongela entero y el daño vuelve a crecer.
    @test count(sol.activos[7]) == 0
    @test maximum(sol.alpha[7]) > maximum(sol.alpha[3]) + 1e-3
    # La energía reducida no crece en ningún barrido de la alternada.
    r = alternancia_at2(zeros(p.N + 1), p, 0.9 * εc * p.L)
    # Dos barridos, no más: en la barra homogénea el paso en α da el α uniforme EXACTO de una vez y el
    # segundo barrido solo comprueba que no se movió. Medido el 14-nov; con defecto son decenas.
    @test length(r.energias) == 2
    @test length(alternancia_at2(zeros(p.N + 1), ParAT2(; ℓ = 0.3, N = 60, defecto = 0.3),
                                 0.9 * εc * p.L).energias) > 10
    @test maximum(diff(r.energias)) <= 1e-12 * maximum(abs, r.energias)
    # Control de mutación: sin la cota inferior el paso en α de la descarga SÍ devuelve daño, así que
    # lo que congela α es la irreversibilidad y no que el problema no tenga por dónde bajar.
    A3 = sol.alpha[3]
    H, b = qp_alpha(respuesta_u(A3, p, carga[5]), p)
    libre = copy(A3)
    pgs!(libre, H, b, zeros(p.N + 1), fill(p.amax, p.N + 1))
    @test minimum(libre .- A3) < -1e-4
end

@testset "4.1 AT2: con cotas naturales en α la localización se va al BORDE, no al defecto" begin
    # Hallazgo medido el 14-nov, no predicción: el nodo del borde tiene la mitad de la masa disipativa
    # de un nodo interior y un solo vecino que lo penalice en α', así que es el punto más débil de la
    # barra y se rompe ahí aunque el defecto esté en el centro. Queda como prueba de regresión porque
    # es la razón por la que las pruebas de localización fijan α = 0 en los extremos.
    for (fb, esperado) in ((false, :borde), (true, :centro))
        q = ParAT2(; ℓ = 0.3, N = 60, defecto = 0.3, fija_borde = fb)
        s = resolver_at2(collect(range(0.1, 1.3; length = 25)) .* eps_c(q) .* q.L, q)
        i = argmax(s.alpha[end])
        if esperado == :borde
            @test i == 1 || i == q.N + 1
        else
            @test abs(i - (q.N ÷ 2 + 1)) <= 2
            @test s.kkt < 1e-10
        end
    end
end

@testset "4.1 AT2: la carga de ruptura converge en malla y escala como ℓ^(-1/2)" begin
    pico(q) = (εc = eps_c(q);
               maximum(resolver_at2(collect(range(0.1, 1.3; length = 30)) .* εc .* q.L, q).sigma))
    base = (; L = 1.0, ℓ = 0.25, defecto = 0.4, ancho = 0.1, fija_borde = true)
    ps = [pico(ParAT2(; base..., N = n)) for n in (40, 80, 160)]
    d1, d2 = abs(ps[2] - ps[1]), abs(ps[3] - ps[2])
    # Criterio fijado ANTES de medir: al duplicar la malla el error se divide por al menos 1,8, o sea
    # orden ≥ 0,85. Medido: razón 4,00, orden 2,00 —P1 da orden 2 en la carga de ruptura—, y el orden
    # observado queda en el json del objetivo 5. El criterio no se retocó para calzar con lo medido.
    @test d2 < d1 / 1.8
    # ⚠️ Acá había un supuesto FALSO, y queda escrito porque el error es instructivo: «un pico
    # localizado queda por debajo del de la rama homogénea». No: medido el 14-nov, 0,7709 contra
    # 0,6495. Fijar α = 0 en los extremos es una RESTRICCIÓN extra, así que sube la energía y la barra
    # con bordes reforzados es MÁS resistente que la rama homogénea libre, defecto incluido. Lo que sí
    # se puede exigir —y es lo que dice «localiza»— es que el daño se concentre en el defecto: en el
    # centro es varias veces el de los cuartos.
    q = ParAT2(; base..., N = 160)
    af = resolver_at2(collect(range(0.1, 1.3; length = 30)) .* eps_c(q) .* q.L, q).alpha[end]
    # Perfil medido (N=160): 0 → 0,155 (x=⅛) → 0,315 (x=¼) → 0,546 → 0,999 (centro). La banda es
    # ancha porque ℓ es un cuarto de L, así que la razón contra el cuarto es 3,17 y contra el
    # octavo 6,46; los umbrales son ésos menos un margen.
    @test af[q.N ÷ 2 + 1] > 2.5 * af[q.N ÷ 4 + 1]
    @test af[q.N ÷ 2 + 1] > 5 * af[q.N ÷ 8 + 1]
    @test ps[3] > sigma_pico_homogenea(ParAT2(; ℓ = 0.25))
    # Invariancia exacta de escala: con L/ℓ, h/ℓ y ancho/ℓ fijos, σ_pico·√(ℓ/(E G_c)) no se mueve. Es
    # la ley ℓ^(-1/2) medida sin ajustar ninguna recta.
    a = ParAT2(; L = 1.0, ℓ = 0.25, defecto = 0.4, ancho = 0.1, N = 80, fija_borde = true)
    b = ParAT2(; L = 2.0, ℓ = 0.50, defecto = 0.4, ancho = 0.2, N = 80, fija_borde = true)
    pa, pb = pico(a), pico(b)
    @test abs(pa * sqrt(a.ℓ) - pb * sqrt(b.ℓ)) < 1e-9 * pa * sqrt(a.ℓ)
    # Control de mutación: sin el √ℓ los dos números difieren en un 29 %, así que la invariancia no es
    # la trivialidad de comparar algo consigo mismo.
    @test abs(pa - pb) > 0.2 * pa
end

# ---- F4.2 · el adjunto discreto y el gradiente -------------------------------------------------

@testset "4.2 AT2: el residuo y el Hessiano REDUCIDOS son los que el adjunto necesita" begin
    p = ParAT2(; L = 1.0, ℓ = 0.3, N = 40, defecto = 0.0)
    εc = eps_c(p)
    d = 0.9 * εc * p.L
    # Fuera del estado convergido: ahí el residuo es cero y compararlo con cualquier cosa pasa igual.
    A = min.(resolver_at2([0.3, 0.6, 0.9] .* εc .* p.L, p).alpha[end] .+
             0.07 .* sin.(range(0, 7; length = p.N + 1)) .^ 2, 0.9)
    H, b = qp_alpha(respuesta_u(A, p, d), p)
    g = residuo_at2(A, p, d)
    @test maximum(abs, g .- (H * A .- b)) < 1e-12 * maximum(abs, g)
    # El Hessiano reducido contra diferencias finitas centradas del residuo, columna por columna.
    Hr = hessiano_at2(A, p, d)
    τ = 1e-7
    for j in (1, 7, 21, p.N + 1)
        e = zeros(p.N + 1)
        e[j] = τ
        fd = (residuo_at2(A .+ e, p, d) .- residuo_at2(A .- e, p, d)) ./ (2τ)
        @test maximum(abs, Hr[:, j] .- fd) < 1e-7 * maximum(abs, fd)
    end
    # Y el término de rango uno —el precio de haber eliminado u— no es ruido: mueve el Hessiano.
    @test maximum(abs, Hr .- hessiano_at2(A, p, d; rango1 = false)) > 1e-5 * maximum(abs, Hr)
end

@testset "4.2 AT2: el gradiente adjunto ES la derivada donde el conjunto activo no cambia" begin
    # ÉSTA es la prueba que decide si el objetivo 5 existe: sin gradiente verificado no hay control
    # óptimo en el modelo de EDP, solo una analogía.
    p = ParAT2(; L = 1.0, ℓ = 0.3, N = 40, defecto = 0.0)
    εc = eps_c(p)
    K = 8
    blancos = collect(range(0.06, 0.40; length = K))
    om = 0.01
    carga = collect(range(0.2, 0.95; length = K)) .* εc .* p.L
    r = gradiente_at2(carga, p; objetivo = blancos, omega = om)
    @test all(==(0), r.congelados)
    @test abs(r.J - costo_at2(resolver_at2(carga, p); objetivo = blancos, omega = om).J) < 1e-14
    costo(c) = costo_at2(resolver_at2(c, p); objetivo = blancos, omega = om).J
    τ = 1e-6 * εc * p.L
    for dir in (normalize(ones(K)), normalize(collect(1.0:K)), [k == 4 ? 1.0 : 0.0 for k in 1:K])
        fd = (costo(carga .+ τ .* dir) - costo(carga .- τ .* dir)) / (2τ)
        @test abs(dot(r.grad, dir) - fd) < 1e-5 * max(1.0, abs(fd))   # medido: 1e-10
    end
    # Control de mutación: con el adjunto del PASO EN α —el que saldría de «derivar el solver», sin el
    # término de rango uno— el gradiente sale con el SIGNO CAMBIADO. Es lo que prueba que la
    # eliminación de u está en el adjunto y no solo en el forward.
    rs = gradiente_at2(carga, p; objetivo = blancos, omega = om, rango1 = false)
    dir = normalize(ones(K))
    fd = (costo(carga .+ τ .* dir) - costo(carga .- τ .* dir)) / (2τ)
    @test abs(dot(rs.grad, dir) - fd) > 0.5 * abs(fd)
    @test dot(rs.grad, dir) * fd < 0
end

@testset "4.2 AT2: donde el conjunto activo CAMBIA el adjunto no es la derivada, y eso es el resultado" begin
    # La etapa 4 sostiene la carga de la 3. Ahí α está en la cota de irreversibilidad con multiplicador
    # cero: subir d₄ hace crecer el daño, bajarlo lo deja congelado. La aplicación carga→estado tiene
    # un quiebre y no hay derivada. El adjunto con el conjunto activo congelado da la derivada de la
    # RAMA, que no es la misma cosa, y esta prueba está escrita para medir el desacuerdo en vez de
    # taparlo: es la evidencia numérica de que la aplicación carga→estado no es diferenciable, y por
    # eso el PMP clásico no está disponible. ⚠️ CORREGIDO 14-nov: antes decía «en EDP solo se tiene
    # C-estacionariedad», atribución FALSA; ese nombre es de Hintermüller-Kopacka 2009, no de
    # Neitzel-Wick-Wollner. Ver docs/marco/objetivo5-phase-field.md § 4.2.
    p = ParAT2(; L = 1.0, ℓ = 0.3, N = 40, defecto = 0.0)
    εc = eps_c(p)
    carga = [0.4, 0.6, 0.8, 0.8, 0.9] .* εc .* p.L
    obj, om = 0.25, 0.2
    r = gradiente_at2(carga, p; objetivo = obj, omega = om)
    @test r.congelados[4] > 0            # el conjunto activo de la etapa 4 no es vacío ni es todo
    @test r.congelados[4] < p.N + 1
    @test all(==(0), r.congelados[[1, 2, 3, 5]])
    costo(c) = costo_at2(resolver_at2(c, p); objetivo = obj, omega = om).J
    e4 = [0.0, 0, 0, 1.0, 0]
    base = costo(carga)
    razones = Float64[]
    for f in (1e-5, 1e-7)
        τ = f * εc * p.L
        der = (costo(carga .+ τ .* e4) - base) / τ
        izq = (base - costo(carga .- τ .* e4)) / τ
        cen = (der + izq) / 2
        # Las dos derivadas laterales difieren por un factor grande (medido: 32), y el adjunto no es
        # ninguna de las dos ni su promedio.
        @test izq > 3 * der
        # ⚠️ CORREGIDO 14-nov ~21:30. Antes era `@test abs(r.grad[4] - cen) > 0.05 * abs(cen)`, que
        # se pone ROJO (0,0492 < 0,05) al perturbar ℓ en 3,3e-12 relativo: pinchaba el desempate de
        # redondeo del conjunto activo, no el quiebre. Separación a LAS DOS laterales, que aguanta las
        # cinco corridas de docs/marco/re-medicion-impugnados/ (mínimo 1,507).
        @test min(abs(r.grad[4] / der), abs(izq / r.grad[4])) > 1.25
        push!(razones, izq / der)
    end
    # Y es un QUIEBRE, no ruido numérico: al bajar τ cien veces la razón entre las derivadas laterales
    # no se mueve. Un error de truncamiento sí se movería.
    @test abs(razones[1] - razones[2]) < 0.01 * razones[1]
    # Donde el conjunto activo NO cambia, el mismo código concuerda: el desacuerdo es del punto, no
    # del adjunto. (Control de que esta prueba no está midiendo un adjunto roto.)
    limpia = [0.4, 0.6, 0.8, 0.85, 0.9] .* εc .* p.L
    r2 = gradiente_at2(limpia, p; objetivo = obj, omega = om)
    @test all(==(0), r2.congelados)
    τ = 1e-6 * εc * p.L
    fd = (costo_at2(resolver_at2(limpia .+ τ .* e4, p); objetivo = obj, omega = om).J -
          costo_at2(resolver_at2(limpia .- τ .* e4, p); objetivo = obj, omega = om).J) / (2τ)
    @test abs(r2.grad[4] - fd) < 1e-5 * max(1.0, abs(fd))
end

@testset "4.2 AT2: el control óptimo de seguimiento llega al patrón de forma cerrada" begin
    p = ParAT2(; L = 1.0, ℓ = 0.3, N = 32, defecto = 0.0)
    εc = eps_c(p)
    K = 6
    blancos = collect(range(0.06, 0.38; length = K))
    # El patrón: por rate-independencia, con omega = 0 el óptimo es invertir Λ(d) etapa por etapa.
    dex = carga_exacta_at2(blancos, p)
    sx = resolver_at2(dex, p)
    @test maximum(abs, sx.longitud .- blancos) < 1e-10
    @test costo_at2(sx; objetivo = blancos, omega = 0.0).J < 1e-20
    # Y el descenso proyectado sobre los incrementos lo encuentra sin saberlo.
    o = optimizar_at2(fill(dex[end] / K * 0.6, K), p; objetivo = blancos, omega = 0.0, iter = 300)
    @test all(diff(o.historia) .< 0)                       # el costo baja monótonamente
    @test o.J < 1e-4 * o.historia[1]          # medido: baja por un factor 1,7e4 en 300 pasos
    @test maximum(abs, o.sol.longitud .- blancos) < 0.025 * blancos[end]   # criterio del plan: < 5 %
    @test all(diff(o.carga) .>= 0)                          # la carga sale monótona, como se pidió
    # Con omega > 0 el costo también baja, pero el óptimo se CORRE del patrón: el término de energía
    # no es decoración. Sin esto, «el control llega al blanco» no distingue entre los dos problemas.
    o2 = optimizar_at2(fill(dex[end] / K * 0.6, K), p; objetivo = blancos, omega = 0.02, iter = 200)
    @test all(diff(o2.historia) .< 0)
    @test maximum(abs, o2.carga .- dex) > 10 * maximum(abs, o.carga .- dex)
end

# ### cifras_informe.jl · las cifras del informe final
# No se revisa cifra por cifra (son cientos): se revisa la COBERTURA (cuántas salieron), que los nombres sean
# usables en LaTeX y únicos, y que una cifra conocida traiga el número del json y no otro.
@testset "cifras del informe: cobertura, nombres y una cifra contra su json" begin
    c, fuera = cifras_informe(; detalle = true)
    @test length(c) > 0
    @test length(c) == 385                      # 374 de los cinco json de final/ + 11 del bloque `borde`
                                                # de avance/resultados/objetivo4_hjb.json (H14), 15-nov
    @test length(fuera) == 71                   # strings, booleanos, nulos y vectores de más de 4 componentes
    nombres = first.(c)
    @test allunique(nombres)
    @test all(n -> occursin(r"^[A-Za-z]+$", n), nombres)   # \newcommand no acepta dígitos ni guiones
    # La regla clave→nombre, con los dos casos que la tuercen: dígitos y separadores.
    @test trozo_nombre("J_directo_N400") == "JDirectoNCuatroCeroCero"
    @test trozo_nombre("reposo_n1") == "ReposoNUno"
    # Una cifra conocida: J del disparo en el acuerdo de los tres métodos, con 4 cifras significativas y coma.
    d = Dict(c)
    b = JSON.parsefile(joinpath(RES_FINAL, "acuerdo_tres_metodos.json"))
    @test round(b["J_disparo"], sigdigits = 4) == 0.3719
    @test d["acuJDisparo"] == "0{,}3719"
    @test d["dosNTramos"] == "50"               # entero, sin decimales
    # H14: el borde artificial de HJB, leído del json del avance. dV0 == 0 exacto es el resultado.
    bh = JSON.parsefile(joinpath(RES_AVANCE, "objetivo4_hjb.json"))["borde"]
    @test bh["dV0"] == 0.0
    @test d["hjbBordeDVCero"] == "0"
    @test bh["lmax"] == 3.015
    @test d["hjbBordeLmax"] == "3{,}015"
    # Y el archivo que se escribe: cabecera de «no editar a mano» y una línea por cifra.
    tmp = joinpath(mktempdir(), "cifras-final.tex")
    escribir_cifras(c, tmp; etiqueta = "prueba")
    lineas = readlines(tmp)
    @test length(lineas) == length(c) + 1
    @test occursin("No editar a mano", lineas[1])
    @test lineas[2] == "\\newcommand{\\cf$(nombres[1])}{$(last(c[1]))}"
end

end;
