# HJB por un esquema semi-Lagrangiano en una malla de (ℓ, P) y pasos de tiempo dt: hacia atrás desde V(T) = Φ,
# V(x, t) = min_U [costo del paso + V(x + dt f, t + dt)], con V interpolado bilinealmente en el pie del paso.

# Celda y coordenada local de X en el eje e. Fuera del eje se extrapola o se recorta al borde.
function ubicar(e, X, modo)
    modo in (:extrapolar, :recortar) || throw(ArgumentError(string(modo)))
    X = modo == :recortar ? clamp(X, e[1], e[end]) : X
    s = (X - e[1]) / (e[2] - e[1])
    i = clamp(floor(Int, s), 0, length(e) - 2)
    i + 1, s - i
end

evaluar(V, (i, s, j, r)) = (1 - s) * (1 - r) * V[i, j] + s * (1 - r) * V[i+1, j] + (1 - s) * r * V[i, j+1] + s * r * V[i+1, j+1]
interpolar(V, ejes, X0, X1, bordes) = evaluar(V, (ubicar(ejes[1], X0, bordes[1])..., ubicar(ejes[2], X1, bordes[2])...))

# Pie de la característica: un paso de punto medio desde x con control u.
function pie(f, x, u, dt)
    k1 = f(x, u)
    k2 = f((x[1] + dt / 2 * k1[1], x[2] + dt / 2 * k1[2]), u)
    x[1] + dt * k2[1], x[2] + dt * k2[2]
end

# Lo que se reutiliza en cada paso de tiempo: dónde cae el pie y el costo del paso por trapecio.
function tabla(f, costo, x, u, dt, ejes, bordes)
    y = pie(f, x, u, dt)
    (ubicar(ejes[1], y[1], bordes[1])..., ubicar(ejes[2], y[2], bordes[2])..., dt / 2 * (costo(x, u) + costo(y, u)))
end

# Como la dinámica es autónoma, los pies se calculan una vez y se reutilizan en todos los pasos; el mínimo se toma
# sobre una lista finita de controles. El where {F, C} obliga a Julia a especializar en f y costo: sin él, el bucle
# queda con despacho dinámico y va unas 25 veces más lento.
function semi_lagrangiano(f::F, costo::C, g, ejes, controles, T, nt, bordes) where {F, C}
    dt = T / nt
    nodos = collect(Iterators.product(ejes...))
    tab = reshape([tabla(f, costo, x, u, dt, ejes, bordes) for u in controles, x in nodos], length(controles), :)
    V = Vector{Matrix{Float64}}(undef, nt + 1)
    pol = Vector{Matrix{Float64}}(undef, nt)
    V[nt+1] = g.(nodos)
    for k in nt-1:-1:0
        Vsig = V[k+2]
        Vk, Pk = similar(Vsig), similar(Vsig)
        # El primer mínimo, como findmin: desigualdad estricta.
        @inbounds for j in axes(tab, 2)
            mejor, a = Inf, 1
            for i in axes(tab, 1)
                e = tab[i, j]
                q = e[5] + evaluar(Vsig, e)
                if q < mejor
                    mejor, a = q, i
                end
            end
            Vk[j] = mejor
            Pk[j] = controles[a]
        end
        V[k+1], pol[k+1] = Vk, Pk
    end
    V, pol
end

function solucion(f, costo, g, ejes, controles, T, nt, bordes)
    V, pol = semi_lagrangiano(f, costo, g, ejes, controles, T, nt, bordes)
    (; V, politica = pol, f, costo, ejes, controles = collect(float.(controles)), T, nt, dt = T / nt, bordes)
end

# Retroalimentación U(t, x): el mismo mínimo del esquema, pero desde x y con el tiempo δ = t_{k+1} - t que falta
# hasta el nodo siguiente. En los nodos δ = dt y coincide con la política guardada.
function feedback(sol, t, x)
    0 <= t < sol.T || throw(ArgumentError("feedback: t = $t fuera de [0, T)"))
    k = min(floor(Int, t / sol.dt + 1e-9), sol.nt - 1)
    δ = (k + 1) * sol.dt - t
    Vsig = sol.V[k+2]
    mejor, umejor = Inf, sol.controles[1]
    for u in sol.controles
        e = tabla(sol.f, sol.costo, x, u, δ, sol.ejes, sol.bordes)
        q = e[5] + evaluar(Vsig, e)
        if q < mejor
            mejor, umejor = q, u
        end
    end
    umejor
end

valor(sol, k, x) = interpolar(sol.V[k+1], sol.ejes, x[1], x[2], sol.bordes)

# En ℓ se extrapola y en P se recorta al borde de la malla.
const BORDES_FRACTURA = (:extrapolar, :recortar)

function hjb_fractura(p, nl, nP, nt, nu; lmax = 2.55, Pmax = 1.2)
    ejes = (range(1.0, lmax, length = nl), range(0.0, Pmax, length = nP))
    solucion((x, u) -> f(x, u, p), (x, u) -> L(x, u, p), x -> Φ(x, p), ejes, range(0.0, p.qmax, length = nu),
             p.T, nt, BORDES_FRACTURA)
end

# Lazo cerrado: el control se recalcula con la ley al inicio de cada uno de n_sub subintervalos y se mantiene
# fijo mientras RK4 integra estado y costo corriente.
function simular_lazo(f, costo, g, ley, x0, T; n_sub = 400, n_rk = 5)
    h = T / n_sub
    dh = h / n_rk
    x = SVector{2, Float64}(x0)
    corriente = 0.0
    ts, Us, Xs = Float64[], Float64[], [x]
    campo(y, u) = push(SVector{2}(f(x2(y), u)), costo(x2(y), u))
    for m in 0:n_sub-1
        t = m * h
        u = ley(t, x)
        y = push(x, 0.0)
        for _ in 1:n_rk
            k1 = campo(y, u)
            k2 = campo(y + dh / 2 * k1, u)
            k3 = campo(y + dh / 2 * k2, u)
            k4 = campo(y + dh * k3, u)
            y = y + dh / 6 * (k1 + 2k2 + 2k3 + k4)
        end
        x = x2(y)
        corriente += y[3]
        push!(ts, t)
        push!(Us, u)
        push!(Xs, x)
    end
    (t = ts, U = Us, X = Xs, J = corriente + g(x), ℓT = x[1], PT = x[2])
end

lazo_cerrado(p, sol; kw...) = simular_lazo((x, u) -> f(x, u, p), (x, u) -> L(x, u, p), x -> Φ(x, p),
                                          (t, x) -> feedback(sol, t, x), [p.ℓ0, p.P0], p.T; kw...)
