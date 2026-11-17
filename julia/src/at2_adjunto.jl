# F4.2 · El adjunto discreto del AT2 y el gradiente respecto del programa de carga.
#
# Es la pieza que decide si el objetivo 5 existe: sin gradiente no hay control óptimo en el modelo de
# EDP, solo una analogía. El forward de `at2.jl` deja, para cada etapa k,
#
#     ∇Ê(α^k; d_k) = 0   en los nodos LIBRES,        α^k_i = α^{k-1}_i   en los CONGELADOS,
#
# con Ê la energía ya reducida en u. Eso es un sistema G^k(α^k, α^{k-1}, d_k) = 0 por etapa, y el
# adjunto se arma con esas mismas filas y se recorre hacia atrás, igual que en el modelo reducido del
# objetivo 3. Dos cosas que no son obvias y que las pruebas miden:
#
#  · El Hessiano reducido ∇²Ê NO es el del paso en α. Eliminar u agrega un término de RANGO UNO
#    (2c·∇S∇Sᵗ/S³) que lo vuelve denso: es el precio de que la tensión sea global en la barra. Con el
#    Hessiano reducido SIN su término de rango uno (rango1=false, que es un control de MUTACIÓN y no
#    el paso en α del solver: ese es qp_alpha, a u fijo, y sus aportes son +k/4) el gradiente no concuerda con
#    diferencias finitas, y eso está puesto como control de mutación.
#  · Congelar el conjunto activo da la derivada de la RAMA, no de la aplicación carga→estado, que en
#    un punto donde el conjunto activo cambia no es derivable. Ese desacuerdo no es un defecto del
#    adjunto: es la evidencia numérica de que la aplicación carga→estado no es diferenciable, y por
#    eso la prueba que lo mide está escrita para que FALLE la comparación y pase la del quiebre.
#    ⚠️ CORREGIDO el 14-nov-2026: antes esta línea decía «en EDP solo se tiene C-estacionariedad»,
#    atribuido a Neitzel-Wick-Wollner. La atribución es FALSA: medido sobre el preprint completo,
#    `C-station`, `M-station` y `weak station` dan 0, con `stationar`=6 y `complementarity`=4 como
#    control positivo. Lo que esa fuente sí dice es que el problema es un MPCC y que las CQ estándar
#    no se satisfacen (p. 1 del preprint); el nombre C-estacionariedad es de su ref. [25],
#    Hintermüller-Kopacka 2009, que NO está en disco. Detalle en docs/marco/objetivo5-phase-field.md.

"∂S/∂A_i con S = Σ_e (1-α_e)^{-2}, que es Σ_{e∋i} (1-α_e)^{-3}."
function grad_S(u, p::ParAT2)
    g = zeros(p.N + 1)
    for e in 1:(p.N)
        v = 1 / u.w[e]^3
        g[e] += v
        g[e + 1] += v
    end
    g
end

"""
    residuo_at2(A, p, d)

∇Ê(A; d) = H_disipado·A − c·∇S/S², con c = E d²/(2h). Es la MISMA cantidad que el gradiente del paso
en α evaluado en la u convergida, escrita de otra forma: la primera es la que usa el adjunto y la
segunda la que usa el solver, y que coincidan es una prueba y no una identidad que se dé por buena.
"""
function residuo_at2(A::AbstractVector, p::ParAT2, d::Real)
    u = respuesta_u(A, p, d)
    c = p.E * d^2 / (2 * paso_malla(p))
    hessiano_disipado(p) * A .- (c / u.S^2) .* grad_S(u, p)
end

"""
    hessiano_at2(A, p, d; rango1)

∇²Ê: la parte tridiagonal (disipada más la elástica local) más el término de rango uno que deja la
eliminación de u. `rango1 = false` lo deja fuera, y existe solo para el control de mutación.
"""
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

"∂²Ê/∂A∂d: la columna por la que la carga entra al sistema adjunto."
dresiduo_dd(A::AbstractVector, p::ParAT2, d::Real) =
    (u = respuesta_u(A, p, d); -(p.E * d / (paso_malla(p) * u.S^2)) .* grad_S(u, p))

"∂(∫α dx)/∂A por trapecios: h en los nodos interiores, h/2 en los extremos."
function grad_longitud(p::ParAT2)
    h = paso_malla(p)
    g = fill(h, p.N + 1)
    g[1] = h / 2
    g[end] = h / 2
    g
end

"""
    trabajo_at2(sol)

Trabajo inyectado ∫σ dd, por trapecios sobre las etapas y con σ(0) = d(0) = 0. Es la energía que el
control gasta, y el término que compite contra llegar al daño pedido.
"""
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

"""
    pesos_objetivo(objetivo, K)

El blanco de daño, por etapa. Un escalar es blanco TERMINAL y un vector es SEGUIMIENTO etapa a etapa.
La distinción no es cosmética y es un resultado del modelo: el AT2 reducido es RATE-INDEPENDENT, así
que con carga monótona α^k es función de d_k y nada más. Con blanco terminal, entonces, J depende solo
de d_K y todos los programas monótonos que lleguen ahí cuestan lo mismo: el problema no determina el
programa. ⚠️ CORREGIDO 14-nov ~21:50: antes decía «con blanco terminal el descenso se estanca en un
valle plano y queda a un 23–33 % del blanco», sin corrida que lo sostuviera —este docstring era su
única fuente—. Medido con `deficit_terminal_at2`, que existe por eso: a omega = 0 el descenso
ALCANZA el blanco terminal (déficit 0,00 % en las dos barras, 12-16 pasos). Lo indeterminado es el
PROGRAMA, no el punto final. El estancamiento aparece con omega > 0 y es del término de trabajo: en
24 combinaciones de barra, omega y punto de partida el déficit fue 0,00 / 0,86 / 5,1 / 7,7-9,2 / 47,8 /
63,2 / 100 %, y ningún 23-33 %. Con seguimiento el programa sí queda determinado, y además tiene
óptimo cerrado para omega = 0 —d_k tal que Λ(d_k) = blanco_k, etapa por etapa—, que es contra lo que
se mide.
"""
pesos_objetivo(objetivo::Real, K::Int) =
    ([k == K ? 1.0 : 0.0 for k in 1:K], fill(Float64(objetivo), K))
pesos_objetivo(objetivo::AbstractVector, K::Int) =
    (length(objetivo) == K ? (ones(K), collect(Float64, objetivo)) :
     throw(ArgumentError("el objetivo de seguimiento necesita una entrada por etapa")))

"""
    costo_at2(sol; objetivo, omega)

J = ½Σ_k peso_k (∫α dx|_k − blanco_k)² + omega·W: «llevar el daño al programa pedido gastando poca
energía», que es el problema de control de la tarea 4.2 y el análogo en daño del LQ de seguimiento de
F1.4. `objetivo` escalar = blanco terminal; vector = seguimiento (ver `pesos_objetivo`).
"""
function costo_at2(sol; objetivo, omega)
    pe, bl = pesos_objetivo(objetivo, length(sol.carga))
    W = trabajo_at2(sol)
    seg = 0.5 * sum(pe[k] * (sol.longitud[k] - bl[k])^2 for k in eachindex(pe))
    (J = seg + omega * W, longitud = sol.longitud[end], trabajo = W, seguimiento = seg)
end

"""
    gradiente_at2(carga, p; objetivo, omega, sol, rango1)

dJ/d(carga) por el adjunto discreto, con el conjunto activo congelado tal como lo dejó el forward.
Devuelve además cuántas filas quedaron congeladas por etapa: con la carga estrictamente creciente son
cero, y ahí el gradiente es la derivada; donde no son cero, es la derivada de la rama.
"""
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
    # Pesos del trabajo: ∂W/∂σ_k y la parte EXPLÍCITA de ∂W/∂d_k (la implícita va por el adjunto).
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

"""
    optimizar_at2(v0, p; objetivo, omega, iter, paso0, tol)

Descenso proyectado sobre los INCREMENTOS v ≥ 0 de la carga, así d = cumsum(v) sale monótona y el
conjunto activo de la irreversibilidad queda vacío, que es justo donde el gradiente adjunto es la
derivada. Paso por retroceso; devuelve la historia del costo, que es lo que se chequea.
"""
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

"""
    carga_exacta_at2(blancos, p; tol)

El óptimo cerrado del seguimiento con omega = 0, para la barra HOMOGÉNEA y blancos crecientes: por
rate-independencia Λ(d) es una función creciente de d, así que basta invertirla etapa por etapa por
bisección. Es el patrón contra el que se mide `optimizar_at2`, no una pieza del método.
"""
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
