# Método directo: el caudal U es constante en N tramos iguales de [0, T] y se optimiza con Ipopt,
# con el gradiente exacto del esquema discreto calculado por el adjunto.
malla(p, N) = range(0.0, p.T, length = N + 1)

# Integra estado y costo corriente tramo a tramo; devuelve J = Φ + ∫L y los tramos con su salida densa.
function avance(q, p)
    t = malla(p, length(q))
    m = cld(PASOS_RK, length(q))
    y = SVector(p.ℓ0, p.P0, 0.0)
    tramos = []
    for k in eachindex(q)
        s = rk4((z, _) -> push(f(z, q[k], p), L(z, q[k], p)), y, (t[k], t[k+1]), m)
        push!(tramos, s)
        y = s.u[end]
    end
    Φ(y, p) + y[3], identity.(tramos)
end

# Integra hacia atrás ζ' = -∂fᵀζ + ∂L desde ζ(T) = -∂Φ, con el mismo RK4. La tercera componente acumula
# -∂H/∂U = wq U - ζP/D en cada tramo, que es la derivada de J respecto de ese tramo del caudal.
function adjunto(q, p, tramos)
    t = malla(p, length(q))
    m = cld(PASOS_RK, length(q))
    xT = tramos[end].u[end]
    ζ = SVector(-p.wT * (xT[1] - p.ℓd), 0.0)
    g = zeros(length(q))
    ζs = Vector{Any}(undef, length(q))
    for k in reverse(eachindex(q))
        xk = tramos[k]
        s = rk4(push(ζ, 0.0), (t[k+1], t[k]), m) do y, τ
            x = xk(τ)
            push(-∂f(x, q[k], p)' * x2(y) + ∂L(x, p), p.wq * q[k] - y[2] / D(x, p))
        end
        g[k] = -s.u[end][3]
        ζ = x2(s.u[end])
        ζs[k] = s
    end
    g, identity.(ζs)
end

costo_y_gradiente(q, p) = (r = avance(q, p); (r[1], adjunto(q, p, r[2])[1]))
gradiente(q, p) = costo_y_gradiente(q, p)[2]

# El problema como NLPModel para Ipopt: el caudal es q = Bθ, con B la identidad o una base (por ejemplo, rampas),
# y cotas 0 ≤ θ ≤ qmax. Ipopt aproxima el Hessiano con L-BFGS.
struct Directo{M} <: AbstractNLPModel{Float64, Vector{Float64}}
    meta::NLPModelMeta{Float64, Vector{Float64}}
    counters::Counters
    p::Par
    B::M
end

Directo(p, B, θ0) = Directo(NLPModelMeta(size(B, 2); x0 = float.(θ0), lvar = zeros(size(B, 2)),
                                         uvar = fill(p.qmax, size(B, 2)), hess_available = false), Counters(), p, B)
NLPModels.obj(m::Directo, θ::AbstractVector) = avance(m.B * θ, m.p)[1]
NLPModels.grad!(m::Directo, θ::AbstractVector, g::AbstractVector) = (g .= m.B' * gradiente(m.B * θ, m.p))

const IPOPT = (print_level = 0, tol = 1e-12, limited_memory_max_history = 30, max_iter = 3000, honor_original_bounds = "yes")

const TOL_ACTIVA = 1e-6

# Ipopt es de punto interior y deja las cotas activas un poco adentro; las que están a menos de TOL_ACTIVA y cuyo
# gradiente empuja hacia afuera se llevan justo a la cota.
function fijar_activas(θ, g, hi)
    arriba = findall(k -> 0 < hi - θ[k] ≤ TOL_ACTIVA && g[k] < 0, eachindex(θ))
    abajo = findall(k -> 0 < θ[k] ≤ TOL_ACTIVA && g[k] > 0, eachindex(θ))
    θ = copy(θ)
    θ[arriba] .= hi
    θ[abajo] .= 0.0
    θ, length(arriba) + length(abajo)
end

function ipopt_directo(p, B, θ0; kw...)
    st = ipopt(Directo(p, B, θ0); IPOPT..., kw...)
    θ, fijadas = fijar_activas(st.solution, B' * gradiente(B * st.solution, p), p.qmax)
    (; x = θ, fun = avance(B * θ, p)[1], fun_sin_fijar = avance(B * st.solution, p)[1], fijadas, exito = st.status == :first_order,
       status = st.status, iter = st.iter)
end

optimizar(p, N; q0 = zeros(N), kw...) = ipopt_directo(p, I(N), q0; kw...)

# Chequeos de optimalidad de un control discreto: norma del gradiente, H casi constante en el tiempo y la brecha
# entre q y U* = proy(ζP / (wq D)) en el punto medio de cada tramo.
function diagnosticos(q, p)
    t = malla(p, length(q))
    J, tramos = avance(q, p)
    g, ζs = adjunto(q, p, tramos)
    m = (t[1:end-1] .+ t[2:end]) ./ 2
    x = [tramos[k](m[k]) for k in eachindex(q)]
    ζ = [ζs[k](m[k]) for k in eachindex(q)]
    Hs = [H(x[k], ζ[k], q[k], p) for k in eachindex(q)]
    brecha = [q[k] - ustar(x[k], ζ[k], p) for k in eachindex(q)]
    Pmin = minimum(minimum(tramos[k](s)[2] for s in range(t[k], t[k+1], length = 9)) for k in eachindex(q))
    xT = tramos[end].u[end]
    Dict("J" => J, "norma_grad_inf" => maximum(abs, g), "H_oscilacion" => maximum(Hs) - minimum(Hs),
         "H_medio" => sum(Hs) / length(Hs), "pmp_brecha_max" => maximum(abs, brecha), "q_min" => minimum(q),
         "q_max" => maximum(q), "q_ultimo" => q[end], "n_cota_inf" => count(≤(1e-12), q),
         "n_cota_sup" => count(≥(p.qmax - 1e-12), q), "P_min" => Pmin, "ell_T" => xT[1], "P_T" => xT[2])
end
