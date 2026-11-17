function integrar(g, y0, tspan; kw...)
    s = solve(ODEProblem{false}((y, _, t) -> g(y, t), y0, tspan), Vern9(); reltol = 1e-12, abstol = 1e-14, kw...)
    s.retcode == ReturnCode.Success || error("integrar: $(s.retcode)")
    s
end

# RK4 de paso fijo con salida densa de Hermite cúbica, solo para avance y adjunto (T5, medido el 13-oct sobre q200):
# PASOS_RK pasos en todo el horizonte, m = cld(PASOS_RK, N) por tramo: en q200 (m = 4), J cambia 5e-14 al duplicarlo
# y difiere 6e-14 del de Vern9. Con m fijo por tramo, el control constante (N = 1) quedaba en 4 pasos. Sirve hacia atrás.
const PASOS_RK = 800

struct TramoRK{V}
    t0::Float64
    h::Float64
    u::Vector{V}
    du::Vector{V}
end

function (s::TramoRK)(τ)
    m = length(s.u) - 1
    i = clamp(floor(Int, (τ - s.t0) / s.h), 0, m - 1)
    θ = (τ - (s.t0 + i * s.h)) / s.h
    h00 = (1 + 2θ) * (1 - θ)^2
    h10 = θ * (1 - θ)^2
    h01 = θ^2 * (3 - 2θ)
    h11 = θ^2 * (θ - 1)
    h00 * s.u[i+1] + (h10 * s.h) * s.du[i+1] + h01 * s.u[i+2] + (h11 * s.h) * s.du[i+2]
end

function rk4(g::G, y0, (ta, tb), m) where {G}
    h = (tb - ta) / m
    y = y0
    u, du = [y0], typeof(y0)[]
    for i in 0:m-1
        t = ta + i * h
        k1 = g(y, t)
        k2 = g(y + h / 2 * k1, t + h / 2)
        k3 = g(y + h / 2 * k2, t + h / 2)
        k4 = g(y + h * k3, t + h)
        y = y + h / 6 * (k1 + 2k2 + 2k3 + k4)
        push!(du, k1)
        push!(u, y)
    end
    push!(du, g(y, tb))
    TramoRK(ta, h, u, du)
end

function simular(p, U, ts::AbstractVector)
    s = integrar((y, _) -> push(f(y, U, p), y[1] * y[2]), SVector(p.ℓ0, p.P0, 0.0), (first(ts), last(ts));
                 saveat = ts)
    ℓ, P, I = (getindex.(s.u, i) for i in 1:3)
    (t = s.t, ℓ, P, S = (p.c .+ ℓ .^ 2) .* P, I)
end

# La malla de salida uniforme es el caso de siempre; la variante de arriba la recibe dada, porque la
# trayectoria libre de la figura del objetivo 1 recorre casi todo su camino en el 1 % inicial del
# horizonte y luego se arrastra al equilibrio: con τ uniforme, el codo del plano de fase queda con
# cuatro puntos y la recta final con setecientos.
simular(p, U, T, npts) = simular(p, U, range(0.0, T, length = npts))
