const N_RAMPA = 200

Jconst(U, p) = avance([U], p)[1]
ℓTconst(U, p) = avance([U], p)[2][end].u[end][1]
info(q, J; extra...) = Dict{String, Any}("J" => J, "U_min" => minimum(q), "U_max" => maximum(q),
                                         (string(k) => v for (k, v) in extra)...)

function constante_factible(p)
    U = raiz(U -> ℓTconst(U, p) - p.ℓd, 0.0, p.qmax)
    U, ℓTconst(U, p)
end

function mejor_constante(p)
    U = raiz(U -> gradiente([U], p)[1], 0.0, p.qmax)
    info([U], Jconst(U, p); U, toca_cota = min(U, p.qmax - U) < 1e-8)
end

function base_rampa(p)
    t = malla(p, N_RAMPA)
    s = (t[1:end-1] .+ t[2:end]) ./ 2 ./ p.T
    [1 .- s s]
end

function mejor_rampa(p)
    B = base_rampa(p)
    U0 = mejor_constante(p)["U"]
    st = ipopt_directo(p, B, [U0, U0])
    a, b = st.x
    info(B * [a, b], st.fun; a, b, escalonada_N = N_RAMPA, toca_cota = any(min.([a, b], p.qmax .- [a, b]) .< 1e-8))
end

function volumen_optimo(p, semilla)
    s = resolver(p, semilla)["sol"]
    quadgk(τ -> (y = s(τ); ustar(x2(y), ζ2(y), p)), 0.0, p.T; rtol = 1e-13)[1]
end

function constante_igual_volumen(p, semilla)
    V = volumen_optimo(p, semilla)
    U = V / p.T
    info([U], Jconst(U, p); U, volumen = V)
end

function desglose(q, p)
    t = malla(p, length(q))
    h = p.T / length(q)
    J, tramos = avance(q, p)
    xT = tramos[end].u[end]
    presion = sum(quadgk(τ -> p.wp * (tramos[k](τ)[2] - p.Pref)^2 / 2, t[k], t[k+1]; atol = 1e-15, rtol = 1e-13)[1]
                  for k in eachindex(q))
    Dict("J" => J, "terminal" => Φ(xT, p), "presion" => presion, "control" => p.wq * sum(abs2, q) * h / 2,
         "volumen" => sum(q) * h, "ell_T" => xT[1], "P_T" => xT[2])
end
