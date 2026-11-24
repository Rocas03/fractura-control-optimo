# Adjunto discreto del AT2 y gradiente del costo respecto del programa de carga (C5).
# En cada etapa k el forward deja ∇Ê(α^k; d_k) = 0 en los nodos libres y α^k_i = α^{k-1}_i en los
# congelados, con Ê la energía ya reducida en u. El adjunto usa esas mismas filas y se recorre hacia atrás,
# igual que en el modelo reducido.
# El Hessiano reducido ∇²Ê no es el del paso en α: eliminar u agrega un término de rango uno,
# 2c·∇S∇Sᵗ/S³, que lo vuelve denso porque la tensión es global en la barra. Sin ese término el gradiente
# no calza con diferencias finitas.
# Congelar el conjunto activo da la derivada de la rama, no la de la aplicación carga → estado, que no es
# derivable donde el conjunto activo cambia. Como señalan Neitzel, Wick y Wollner, el problema es un MPCC y
# las condiciones de calificación estándar no se cumplen.

# ∂S/∂A_i con S = Σ_e (1-α_e)^{-2}, que da Σ_{e∋i} (1-α_e)^{-3}.
function grad_S(u, p::ParAT2)
    g = zeros(p.N + 1)
    for e in 1:(p.N)
        v = 1 / u.w[e]^3
        g[e] += v
        g[e + 1] += v
    end
    g
end

# ∇Ê(A; d) = H_disipado·A − c·∇S/S², con c = E d²/(2h). Es el gradiente del paso en α evaluado en la u
# óptima, escrito de otra forma; las pruebas comprueban que las dos formas coinciden.
function residuo_at2(A::AbstractVector, p::ParAT2, d::Real)
    u = respuesta_u(A, p, d)
    c = p.E * d^2 / (2 * paso_malla(p))
    hessiano_disipado(p) * A .- (c / u.S^2) .* grad_S(u, p)
end

# ∇²Ê: la parte tridiagonal (disipada más elástica local) más el término de rango uno que deja eliminar u.
# `rango1 = false` lo omite y sirve solo como control de mutación en las pruebas.
function hessiano_at2(A::AbstractVector, p::ParAT2, d::Real; rango1 = true)
    u = respuesta_u(A, p, d)
    c = p.E * d^2 / (2 * paso_malla(p))
    H = Matrix(hessiano_disipado(p))
    for e in 1:(p.N)
        v = -c * 1.5 / (u.w[e]^4 * u.S^2)
        H[e, e] += v
        H[e, e + 1] += v
        H[e + 1, e] += v
        H[e + 1, e + 1] += v
    end
    if rango1
        Sg = grad_S(u, p)
        H .+= (2 * c / u.S^3) .* (Sg * Sg')
    end
    H
end

# ∂²Ê/∂A∂d: la columna por la que la carga entra al sistema adjunto.
dresiduo_dd(A::AbstractVector, p::ParAT2, d::Real) =
    (u = respuesta_u(A, p, d); -(p.E * d / (paso_malla(p) * u.S^2)) .* grad_S(u, p))

# ∂(∫α dx)/∂A por trapecios: h en los nodos interiores y h/2 en los extremos.
function grad_longitud(p::ParAT2)
    h = paso_malla(p)
    g = fill(h, p.N + 1)
    g[1] = h / 2
    g[end] = h / 2
    g
end

# Trabajo inyectado W = ∫σ dd, por trapecios sobre las etapas con σ(0) = d(0) = 0: la energía que gasta
# el control, que compite con llegar al daño pedido.
function trabajo_at2(sol)
    d, σ = sol.carga, sol.sigma
    W = 0.0
    dp = 0.0
    σp = 0.0
    for k in eachindex(d)
        W += 0.5 * (σ[k] + σp) * (d[k] - dp)
        σp = σ[k]
        dp = d[k]
    end
    W
end

# Blanco de daño por etapa: un escalar es blanco terminal y un vector es seguimiento etapa a etapa.
# El AT2 reducido no depende de la velocidad de carga, así que con carga monótona α^k depende solo de d_k.
# Con blanco terminal J depende solo de d_K: el descenso llega al blanco, pero el programa queda
# indeterminado. Con seguimiento sí queda determinado, y con omega = 0 tiene óptimo cerrado, Λ(d_k) =
# blanco_k en cada etapa con Λ el largo dañado, que es contra lo que se compara.
pesos_objetivo(objetivo::Real, K::Int) =
    ([k == K ? 1.0 : 0.0 for k in 1:K], fill(Float64(objetivo), K))
pesos_objetivo(objetivo::AbstractVector, K::Int) =
    (length(objetivo) == K ? (ones(K), collect(Float64, objetivo)) :
     throw(ArgumentError("el objetivo de seguimiento necesita una entrada por etapa")))

# J = ½Σ_k peso_k (∫α dx|_k − blanco_k)² + omega·W: llevar el daño al programa pedido gastando poca
# energía, el análogo en daño del LQ de seguimiento. `objetivo` escalar es blanco terminal; vector,
# seguimiento.
function costo_at2(sol; objetivo, omega)
    pe, bl = pesos_objetivo(objetivo, length(sol.carga))
    W = trabajo_at2(sol)
    seg = 0.5 * sum(pe[k] * (sol.longitud[k] - bl[k])^2 for k in eachindex(pe))
    (J = seg + omega * W, longitud = sol.longitud[end], trabajo = W, seguimiento = seg)
end

# dJ/d(carga) por el adjunto discreto, con el conjunto activo congelado tal como lo dejó el forward.
# También devuelve cuántos nodos quedaron congelados por etapa: con carga estrictamente creciente son cero
# y el gradiente es la derivada; si no, es la derivada de la rama.
function gradiente_at2(carga::AbstractVector, p::ParAT2; objetivo, omega, sol = nothing,
                       rango1 = true)
    s = sol === nothing ? resolver_at2(carga, p) : sol
    K = length(s.carga)
    n = p.N + 1
    h = paso_malla(p)
    d, σ = s.carga, s.sigma
    us = [respuesta_u(s.alpha[k], p, d[k]) for k in 1:K]
    dσ_dA = [-(p.E * d[k] / (h * us[k].S^2)) .* grad_S(us[k], p) for k in 1:K]
    dσ_dd = [p.E / (h * us[k].S) for k in 1:K]
    # Pesos del trabajo: ∂W/∂σ_k y la parte explícita de ∂W/∂d_k (la implícita pasa por el adjunto).
    wσ = [0.5 * (d[k] - (k == 1 ? 0.0 : d[k - 1])) + (k < K ? 0.5 * (d[k + 1] - d[k]) : 0.0)
          for k in 1:K]
    wd = [0.5 * (σ[k] + (k == 1 ? 0.0 : σ[k - 1])) - (k < K ? 0.5 * (σ[k + 1] + σ[k]) : 0.0)
          for k in 1:K]
    congelado = [s.activos[k] .| s.tope[k] .| s.fijos for k in 1:K]
    gΛ = grad_longitud(p)
    pe, bl = pesos_objetivo(objetivo, K)
    grad = zeros(K)
    acc = zeros(n)
    for k in K:-1:1
        dJdA = omega * wσ[k] .* dσ_dA[k]
        pe[k] == 0 || (dJdA .+= (pe[k] * (s.longitud[k] - bl[k])) .* gΛ)
        M = hessiano_at2(s.alpha[k], p, d[k]; rango1)
        for i in 1:n
            if congelado[k][i]
                M[i, :] .= 0.0
                M[i, i] = 1.0
            end
        end
        λ = transpose(M) \ (-(dJdA .+ acc))
        rd = dresiduo_dd(s.alpha[k], p, d[k])
        aporte = 0.0
        for i in 1:n
            congelado[k][i] || (aporte += λ[i] * rd[i])
        end
        grad[k] = omega * (wd[k] + wσ[k] * dσ_dd[k]) + aporte
        acc .= 0.0
        for i in 1:n
            s.activos[k][i] && (acc[i] = -λ[i])
        end
    end
    (J = costo_at2(s; objetivo, omega).J, grad = grad, sol = s,
     congelados = [count(c) for c in congelado])
end

# Descenso proyectado sobre los incrementos v ≥ 0 de la carga: d = cumsum(v) sale monótona y el conjunto
# activo de la irreversibilidad queda vacío, que es donde el gradiente adjunto es la derivada. El paso se
# ajusta por retroceso; devuelve la historia del costo.
function optimizar_at2(v0::AbstractVector, p::ParAT2; objetivo, omega, iter = 60, paso0 = 1.0,
                       tol = 1e-12)
    v = max.(copy(v0), 0.0)
    r = gradiente_at2(cumsum(v), p; objetivo, omega)
    hist = [r.J]
    paso = paso0
    for _ in 1:iter
        g = [sum(@view r.grad[j:end]) for j in eachindex(v)]
        bajo = false
        for _ in 1:40
            w = max.(v .- paso .* g, 0.0)
            rw = gradiente_at2(cumsum(w), p; objetivo, omega)
            if rw.J < r.J
                v, r = w, rw
                push!(hist, r.J)
                paso *= 2
                bajo = true
                break
            end
            paso /= 2
        end
        bajo || break
        abs(hist[end - 1] - hist[end]) < tol * max(1.0, abs(hist[end])) && break
    end
    (v = v, carga = cumsum(v), J = r.J, historia = hist, sol = r.sol, grad = r.grad,
     longitud = r.sol.longitud[end], trabajo = trabajo_at2(r.sol))
end

# Óptimo cerrado del seguimiento con omega = 0, para la barra homogénea y blancos crecientes: Λ(d) crece
# con d, así que basta invertirla etapa por etapa por bisección. Es la referencia contra la que se mide
# `optimizar_at2`.
function carga_exacta_at2(blancos::AbstractVector, p::ParAT2; tol = 1e-12)
    Λ(d) = longitud_danada(resolver_at2([d], p).alpha[end], p)
    [begin
         a, b = 0.0, 8 * eps_c(p) * p.L
         while b - a > tol * b
             m = (a + b) / 2
             Λ(m) < blanco ? (a = m) : (b = m)
         end
         (a + b) / 2
     end for blanco in blancos]
end
