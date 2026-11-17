const TOL_RESIDUO = 1e-10

function extremal(z0, p)
    integrar(SVector(p.ℓ0, p.P0, z0[1], z0[2], 0.0), (0.0, p.T)) do y, _
        x, ζ = x2(y), ζ2(y)
        u = ustar(x, ζ, p)
        vcat(f(x, u, p), -∂f(x, u, p)' * ζ + ∂L(x, p), SVector(L(x, u, p)))
    end
end

function residuo(z0, p, variante)
    yT = extremal(z0, p).u[end]
    variante == :penalizado ? SVector(yT[3] + p.wT * (yT[1] - p.ℓd), yT[4]) : SVector(yT[1] - p.ℓd, yT[4])
end

function conmutaciones(s, p; npts = 4001)
    t = range(0.0, p.T, length = npts)
    w(τ) = (y = s(τ); ulibre(x2(y), ζ2(y), p))
    W = w.(t)
    out = Float64[]
    for nivel in (0.0, p.qmax), i in 1:npts-1
        sign(W[i] - nivel) * sign(W[i+1] - nivel) < 0 && push!(out, raiz(τ -> w(τ) - nivel, t[i], t[i+1]))
    end
    sort(filter(τ -> p.T - τ > 1e-6, out))
end

function resolver(p, semilla, variante = :penalizado)
    variante in (:penalizado, :terminal) || throw(ArgumentError(string(variante)))
    r = newton2(z -> residuo(z, p, variante), semilla)
    s = extremal(r.u, p)
    yT = s.u[end]
    J = yT[5] + (variante == :penalizado ? Φ(yT, p) : 0.0)
    Hs = [(y = s(τ); H(x2(y), ζ2(y), ustar(x2(y), ζ2(y), p), p)) for τ in range(0.0, p.T, length = 401)]
    res = residuo(r.u, p, variante)
    Dict("exito" => maximum(abs, res) < TOL_RESIDUO, "mensaje" => r.retcode, "z0" => r.u, "J" => J,
         "ell_T" => yT[1], "P_T" => yT[2], "zeta_ell_T" => yT[3], "norma_residuo" => maximum(abs, res),
         "H_oscilacion" => maximum(Hs) - minimum(Hs), "conmutaciones" => conmutaciones(s, p), "sol" => s)
end
