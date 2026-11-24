# Parámetros adimensionales del caso base: exponente n de la ley de propagación, c, filtración η, estado inicial,
# horizonte T, largo objetivo ℓd, pesos del costo J y cota qmax del caudal U.
Base.@kwdef struct Par
    n::Float64 = 2.0
    c::Float64 = 1.0
    η::Float64 = 0.2
    ℓ0::Float64 = 1.0
    P0::Float64 = 0.3
    T::Float64 = 2.0
    ℓd::Float64 = 1.7
    wT::Float64 = 10.0
    wp::Float64 = 1.0
    wq::Float64 = 0.2
    qmax::Float64 = 2.0
    Pref::Float64 = 0.0
end

# Copia de p con algunos campos cambiados, y su versión en Dict para escribir los json.
con(p::Par; kw...) = Par(; (k => getfield(p, k) for k in fieldnames(Par))..., kw...)
comodict(p::Par) = Dict("n" => p.n, "c" => p.c, "eta" => p.η, "ell0" => p.ℓ0, "P0" => p.P0, "T" => p.T,
                        "Ld" => p.ℓd, "wT" => p.wT, "wp" => p.wp, "wq" => p.wq, "qmax" => p.qmax, "Pref" => p.Pref)

# D = c + ℓ² es el almacenamiento; x2 y ζ2 separan estado y adjunto en el vector extendido.
D(x, p) = p.c + x[1]^2
x2(y) = SVector(y[1], y[2])
ζ2(y) = SVector(y[3], y[4])

# Dinámica: ℓ' = F = ℓ^(n/2) P^n y D P' = U - ηℓP - 2ℓPF (inyección menos filtración y crecimiento).
function f(x, u, p)
    ℓ, P = x[1], x[2]
    F = ℓ^(p.n / 2) * P^p.n
    SVector(F, (u - p.η * ℓ * P - 2ℓ * P * F) / D(x, p))
end

# Jacobiano de f respecto de (ℓ, P). Con P = 0 se usa n·ℓ^(n/2)·P^(n-1), no nF/P.
function ∂f(x, u, p)
    ℓ, P, n = x[1], x[2], p.n
    d = D(x, p)
    F = ℓ^(n / 2) * P^n
    Fℓ = n / 2 * ℓ^(n / 2 - 1) * P^n
    FP = n * ℓ^(n / 2) * P^(n - 1)
    h = (u - p.η * ℓ * P - 2ℓ * P * F) / d
    Rℓ = p.η * P + 2P * F + 2ℓ * P * Fℓ
    RP = p.η * ℓ + 2ℓ * F + 2ℓ * P * FP
    SMatrix{2, 2}(Fℓ, -Rℓ / d - 2ℓ * h / d, FP, -RP / d)
end

# Costo corriente L, costo terminal Φ y pre-Hamiltoniano H = ζ·f - L, con ζ el adjunto (precio sombra).
# Maximizar H en U da U = ζP / (wq D), proyectado a [0, qmax].
L(x, u, p) = (p.wp * (x[2] - p.Pref)^2 + p.wq * u^2) / 2
∂L(x, p) = SVector(zero(x[2]), p.wp * (x[2] - p.Pref))
Φ(x, p) = p.wT * (x[1] - p.ℓd)^2 / 2
H(x, ζ, u, p) = (v = f(x, u, p); ζ[1] * v[1] + ζ[2] * v[2] - L(x, u, p))
ulibre(x, ζ, p) = ζ[2] / (p.wq * D(x, p))
ustar(x, ζ, p) = clamp(ulibre(x, ζ, p), 0.0, p.qmax)
