const NIVELES = ((41, 41, 10, 21), (81, 81, 20, 41), (161, 161, 40, 81))
const BORDE_AMPLIO = (nl = 105, nP = 101, lmax = 1.0 + 104 * 1.55 / 80, Pmax = 1.5)

function correr_nivel(p, nl, nP, nt, nu; kw...)
    t_hjb = @elapsed sol = hjb_fractura(p, nl, nP, nt, nu; kw...)
    t_lc = @elapsed lc = lazo_cerrado(p, sol)
    sol, lc, t_hjb, t_lc
end

function ∇V(sol, k, x)
    h = [e[2] - e[1] for e in sol.ejes]
    [(valor(sol, k, x .+ h[i] .* (1:2 .== i)) - valor(sol, k, x .- h[i] .* (1:2 .== i))) / (2h[i]) for i in 1:2]
end

function objetivo4()
    p = Par()
    x0 = [p.ℓ0, p.P0]
    rd = resolver(p, SEMILLA)
    Jstar = rd["J"]
    z0 = rd["z0"]
    res = Dict{String, Any}("params" => comodict(p), "J_disparo" => Jstar, "z0_disparo" => z0, "disparo_exito" => rd["exito"],
                            "disparo_norma_residuo" => rd["norma_residuo"], "niveles" => Dict{String, Any}[])
    criterios = Dict{String, Any}("disparo_converge" => rd["exito"])
    sols, lcs, tiempos = [], [], []
    for (nl, nP, nt, nu) in NIVELES
        sol, lc, t_hjb, t_lc = correr_nivel(p, nl, nP, nt, nu)
        i0 = argmin(abs.(sol.ejes[1] .- x0[1]))
        j0 = argmin(abs.(sol.ejes[2] .- x0[2]))
        es_nodo = abs(sol.ejes[1][i0] - x0[1]) < 1e-12 && abs(sol.ejes[2][j0] - x0[2]) < 1e-12
        V0 = sol.V[1][i0, j0]
        Y = rd["sol"].(lc.t)
        us = [ustar(x2(y), ζ2(y), p) for y in Y]
        dz, rel = Vector{Float64}[], Vector{Float64}[]
        for k in 1:nt-1
            yk = rd["sol"](k * sol.dt)
            d = abs.(-∇V(sol, k, x2(yk)) .- ζ2(yk))
            push!(dz, d)
            push!(rel, d ./ abs.(ζ2(yk)))
        end
        dzmax = [maximum(getindex.(dz, i)) for i in 1:2]
        fila = Dict{String, Any}("nl" => nl, "nP" => nP, "nt" => nt, "n_controles" => nu,
                                 "dominio_ell" => [sol.ejes[1][1], sol.ejes[1][end]], "dominio_P" => [sol.ejes[2][1], sol.ejes[2][end]],
                                 "x0_es_nodo" => es_nodo, "V0" => V0, "e_V" => abs(V0 - Jstar), "e_V_rel" => abs(V0 - Jstar) / Jstar,
                                 "J_lc" => lc.J, "dJ_lc" => lc.J - Jstar, "sup_dU" => maximum(abs, lc.U .- us),
                                 "max_U_lc" => maximum(lc.U), "ell_T_lc" => lc.ℓT, "P_T_lc" => lc.PT, "dzeta_max" => dzmax,
                                 "dzeta_max_sobre_z0" => dzmax ./ abs.(z0),
                                 "dzeta_rel_puntual_max" => [maximum(getindex.(rel, i)) for i in 1:2])
        push!(res["niveles"], fila)
        push!(sols, sol)
        push!(lcs, lc)
        push!(tiempos, (nl, t_hjb, t_lc))
    end
    F = res["niveles"]
    e = [f["e_V"] for f in F]
    criterios["x0_es_nodo"] = all(f["x0_es_nodo"] for f in F)
    res["cocientes_e_V"] = [e[i] / e[i+1] for i in 1:length(e)-1]
    criterios["valor_cociente_por_nivel_ge_1.5"] = all(≥(1.5), res["cocientes_e_V"])
    criterios["valor_malla_fina_bajo_1pct"] = e[end] < 0.01 * Jstar
    dU = [f["sup_dU"] for f in F]
    dJ = [abs(f["dJ_lc"]) for f in F]
    criterios["lazo_sup_dU_baja"] = all(dU[i+1] < dU[i] for i in 1:length(dU)-1)
    criterios["lazo_dJ_baja"] = all(dJ[i+1] < dJ[i] for i in 1:length(dJ)-1)
    criterios["lazo_J_no_baja_de_Jstar"] = all(f["J_lc"] ≥ Jstar - 1e-6 for f in F)
    nl, nP, nt, nu = NIVELES[2]
    b = BORDE_AMPLIO
    sol_b, lc_b, t_b, _ = correr_nivel(p, b.nl, b.nP, nt, nu; lmax = b.lmax, Pmax = b.Pmax)
    h_med = [e_[2] - e_[1] for e_ in sols[2].ejes]
    h_b = [e_[2] - e_[1] for e_ in sol_b.ejes]
    V0_b = valor(sol_b, 0, x0)
    res["borde"] = Dict{String, Any}("nl" => b.nl, "nP" => b.nP, "lmax" => b.lmax, "Pmax" => b.Pmax, "h_medio" => h_med,
                                     "h_amplio" => h_b, "V0_amplio" => V0_b, "dV0" => V0_b - F[2]["V0"], "J_lc_amplio" => lc_b.J)
    push!(tiempos, (b.nl, t_b, 0.0))
    criterios["borde_mismo_paso"] = all(abs.(h_med .- h_b) .≤ 1e-12)
    criterios["borde_dV0_bajo_1e-4"] = abs(V0_b - F[2]["V0"]) < 1e-4
    criterios["lazo_U_bajo_qmax"] = all(f["max_U_lc"] < p.qmax for f in F)
    res["criterios"] = criterios
    res, (; rd, lcs, tiempos)
end
