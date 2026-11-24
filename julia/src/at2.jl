# Modelo de daño por gradiente AT2 en una barra 1D (C5), el modelo de EDP al que se traslada el control.
# La energía es E(u, α) = ∫ ½ (1-α)² E u'² dx + (G_c/2) ∫ (α²/ℓ + ℓ α'²) dx, con w(α) = α², sin umbral
# y con la normalización de Pham–Marigo (c_w = 2).
# En 1D sin carga de volumen la tensión es constante en la barra, así que el paso en u sale en forma cerrada.
# AT2 es cuadrático en α, así que el paso en α es un problema cuadrático convexo con cotas; la
# irreversibilidad entra como cota inferior, α mayor o igual que en la etapa de carga anterior.
# Como control sirve la rama homogénea: con cotas naturales, α uniforme es punto estacionario exacto del
# problema discreto, con α(ε) = ε²/(ε² + ε_c²), σ(ε) = E ε (1-α)² y ε_c² = G_c/(E ℓ).
# Para que el daño se localice hay que romper la simetría: `defecto` baja G_c con una gaussiana de ancho
# fijo, que no se encoge al refinar la malla.

struct ParAT2
    L::Float64        # largo de la barra
    E::Float64        # módulo
    Gc::Float64       # tenacidad
    ℓ::Float64        # longitud de regularización
    N::Int            # elementos (P1 en u y en α)
    defecto::Float64  # amplitud del hoyo de G_c, 0 = barra homogénea
    ancho::Float64    # ancho del hoyo, fijo en x: no se encoge al refinar la malla
    amax::Float64     # cota superior de α; si se activa, el modelo degeneró y hay que avisar
    fija_borde::Bool  # α = 0 en los dos nodos extremos
end

# Sobre `fija_borde`: con cotas naturales el nodo del borde es el más débil de la barra (la mitad de masa
# disipativa y un solo vecino), y el daño se localiza en el extremo en vez de en el defecto. Fijar α = 0 en
# los extremos es la condición estándar y la usan las pruebas de localización. La rama homogénea, en
# cambio, necesita cotas naturales: con α fijo en el borde no hay solución uniforme.

ParAT2(; L = 1.0, E = 1.0, Gc = 1.0, ℓ = 0.5, N = 60, defecto = 0.0, ancho = 0.15,
        amax = 1 - 1e-6, fija_borde = false) =
    ParAT2(L, E, Gc, ℓ, N, defecto, ancho, amax, fija_borde)

# Nodos con α fijo en 0: ninguno, o los dos extremos si `fija_borde`.
function nodos_fijos(p::ParAT2)
    f = falses(p.N + 1)
    if p.fija_borde
        f[1] = true
        f[end] = true
    end
    f
end

paso_malla(p::ParAT2) = p.L / p.N
nodos_at2(p::ParAT2) = collect(range(0.0, p.L; length = p.N + 1))
centros_at2(p::ParAT2) = (x = nodos_at2(p); (x[1:(end - 1)] .+ x[2:end]) ./ 2)

# G_c por elemento; con `defecto = 0` es constante y la barra admite la solución homogénea exacta.
gc_elem(p::ParAT2) = p.defecto == 0 ? fill(p.Gc, p.N) :
    [p.Gc * (1 - p.defecto * exp(-((x - p.L / 2) / p.ancho)^2)) for x in centros_at2(p)]

# α en el punto medio de cada elemento, donde se evalúa la degradación elástica.
alpha_elem(A::AbstractVector) = (A[1:(end - 1)] .+ A[2:end]) ./ 2

# ∫α dx por trapecios: cuánto largo de barra está dañado, que es lo que mide el costo del control.
longitud_danada(A::AbstractVector, p::ParAT2) = paso_malla(p) * sum(alpha_elem(A))

# Minimiza en u a α fijo, con u(0) = 0 y u(L) = d. Devuelve σ (constante en la barra), los saltos de u por
# elemento y S = Σ_e (1-α_e)^{-2}, junto con las potencias de (1-α_e) que usan el gradiente y el Hessiano.
function respuesta_u(A::AbstractVector, p::ParAT2, d::Real)
    αe = alpha_elem(A)
    w = 1 .- αe
    s = @. 1 / w^2
    S = sum(s)
    h = paso_malla(p)
    σ = p.E * d / (h * S)
    (σ = σ, S = S, s = s, αe = αe, w = w, du = d .* s ./ S, elastica = 0.5 * σ * d)
end

# (G_c/2)∫(α²/ℓ + ℓα'²), con ∫α² exacta sobre P1 (matriz de masa) y α' constante por elemento.
function energia_disipada(A::AbstractVector, p::ParAT2)
    h = paso_malla(p)
    g = gc_elem(p)
    tot = 0.0
    for e in 1:(p.N)
        a, b = A[e], A[e + 1]
        tot += (g[e] / 2) * ((h / (3 * p.ℓ)) * (a^2 + a * b + b^2) + (p.ℓ / h) * (b - a)^2)
    end
    tot
end

# Hessiano de `energia_disipada`: masa (h/3ℓ, h/6ℓ) más rigidez (ℓ/h, -ℓ/h), ensamblado por elemento.
function hessiano_disipado(p::ParAT2)
    h = paso_malla(p)
    g = gc_elem(p)
    dv = zeros(p.N + 1)
    ev = zeros(p.N)
    for e in 1:(p.N)
        m = g[e] * h / (3 * p.ℓ)
        k = g[e] * p.ℓ / h
        dv[e] += m + k
        dv[e + 1] += m + k
        ev[e] += m / 2 - k
    end
    SymTridiagonal(dv, ev)
end

# El paso en α a u fijo: min ½AᵗHA − bᵗA con α_prev ≤ A ≤ amax. H suma la parte disipada y la elástica,
# y b sale del término lineal de (1-α)². `u` es la salida de `respuesta_u`.
function qp_alpha(u, p::ParAT2)
    h = paso_malla(p)
    Hd = hessiano_disipado(p)
    dv = copy(Hd.dv)
    ev = copy(Hd.ev)
    b = zeros(p.N + 1)
    for e in 1:(p.N)
        k = p.E * h * (u.du[e] / h)^2        # k_e = E h ε_e²
        dv[e] += k / 4
        dv[e + 1] += k / 4
        ev[e] += k / 4
        b[e] += k / 2
        b[e + 1] += k / 2
    end
    (SymTridiagonal(dv, ev), b)
end

# Gauss–Seidel proyectado para ½AᵗHA − bᵗA con cotas: minimiza coordenada por coordenada, así que la
# proyección es un `clamp`. Con H definida positiva cada barrido baja la energía.
function pgs!(A::AbstractVector, H::SymTridiagonal, b::AbstractVector, lo::AbstractVector,
              hi::AbstractVector; iter = 100_000, tol = 1e-15)
    dv, ev = H.dv, H.ev
    n = length(A)
    usadas = iter
    for it in 1:iter
        δ = 0.0
        for i in 1:n
            r = b[i]
            i > 1 && (r -= ev[i - 1] * A[i - 1])
            i < n && (r -= ev[i] * A[i + 1])
            nuevo = clamp(r / dv[i], lo[i], hi[i])
            δ = max(δ, abs(nuevo - A[i]))
            A[i] = nuevo
        end
        if δ < tol
            usadas = it
            break
        end
    end
    (A = A, iteraciones = usadas)
end

# El mismo problema por conjunto activo primal: cada vuelta resuelve el sistema tridiagonal en los nodos
# libres, manda a la cota los que la violan y suelta los que tienen el multiplicador con el signo
# equivocado. En mallas finas es mucho más rápido que Gauss–Seidel. Cierra con un pulido de `pgs!` y
# devuelve el residuo KKT.
function qp_cotas!(A::AbstractVector, H::SymTridiagonal, b::AbstractVector, lo::AbstractVector,
                   hi::AbstractVector; iter = 300)
    n = length(A)
    dv, ev = H.dv, H.ev
    act = falses(n)
    @inbounds for i in 1:n
        A[i] = clamp(A[i], lo[i], hi[i])
    end
    vueltas = iter
    for it in 1:iter
        libres = findall(!, act)
        if !isempty(libres)
            m = length(libres)
            fdv = [dv[i] for i in libres]
            fev = [libres[k + 1] == libres[k] + 1 ? ev[libres[k]] : 0.0 for k in 1:(m - 1)]
            rhs = similar(fdv)
            for k in 1:m
                i = libres[k]
                v = b[i]
                (i > 1 && act[i - 1]) && (v -= ev[i - 1] * A[i - 1])
                (i < n && act[i + 1]) && (v -= ev[i] * A[i + 1])
                rhs[k] = v
            end
            x = SymTridiagonal(fdv, fev) \ rhs
            for k in 1:m
                A[libres[k]] = x[k]
            end
        end
        fuera = false
        for i in 1:n
            act[i] && continue
            if A[i] < lo[i]
                A[i] = lo[i]; act[i] = true; fuera = true
            elseif A[i] > hi[i]
                A[i] = hi[i]; act[i] = true; fuera = true
            end
        end
        fuera && continue
        g = H * A .- b
        # Umbral para soltar un nodo: si la etapa repite la carga anterior, los nodos quedan en la cota con
        # multiplicador cero, y soltarlos por un -1e-16 hace ciclar el conjunto activo.
        umb = 1e-12 * max(1.0, maximum(abs, b))
        suelto = false
        for i in 1:n
            act[i] || continue
            if A[i] <= lo[i] && g[i] < -umb
                act[i] = false; suelto = true
            elseif A[i] >= hi[i] && g[i] > umb
                act[i] = false; suelto = true
            end
        end
        if !suelto
            vueltas = it
            break
        end
    end
    # Si el conjunto activo agotó sus vueltas, salió ciclando y A no cumple KKT; solo en ese caso el pulido
    # pasa a ser el solver, con un tope de barridos mucho mayor.
    pgs!(A, H, b, lo, hi; iter = vueltas == iter ? 200_000 : 50, tol = 1e-15)
    g = H * A .- b
    kkt = 0.0
    for i in 1:n
        # Un nodo con lo == hi está fijo y no tiene condición de optimalidad: el residuo cuenta solo los nodos
        # que el problema puede mover.
        lo[i] >= hi[i] && continue
        if A[i] > lo[i] && A[i] < hi[i]
            kkt = max(kkt, abs(g[i]))
        elseif A[i] <= lo[i]
            kkt = max(kkt, max(0.0, -g[i]))
        else
            kkt = max(kkt, max(0.0, g[i]))
        end
    end
    (A = A, vueltas = vueltas, kkt = kkt)
end

# Minimización alternada en u y α para una etapa de carga, con la irreversibilidad como cota A ≥ A_prev.
# Devuelve α, los barridos y las energías reducidas Ê(A^m) = min_u E(u, A^m), que no deben crecer.
function alternancia_at2(A_prev::AbstractVector, p::ParAT2, d::Real; tol = 1e-13, iter = 2_000)
    A = copy(A_prev)
    hi = fill(p.amax, length(A))
    for i in findall(nodos_fijos(p))
        hi[i] = 0.0
    end
    energias = Float64[]
    usadas = iter
    kkt = 0.0
    for it in 1:iter
        u = respuesta_u(A, p, d)
        H, b = qp_alpha(u, p)
        previo = copy(A)
        q = qp_cotas!(A, H, b, A_prev, hi)
        kkt = max(kkt, q.kkt)
        push!(energias, respuesta_u(A, p, d).elastica + energia_disipada(A, p))
        if maximum(abs, A .- previo) < tol
            usadas = it
            break
        end
    end
    (alpha = A, iteraciones = usadas, energias = energias, kkt = kkt)
end

# Resuelve toda la historia de carga (desplazamiento impuesto en x = L), arrastrando la irreversibilidad
# entre etapas. Para el adjunto guarda qué nodos quedaron pegados a la cota de irreversibilidad (`activos`)
# y a amax (`tope`, que debería quedar vacío).
function resolver_at2(carga::AbstractVector, p::ParAT2 = ParAT2(); tol = 1e-13, iter = 2_000)
    A = zeros(p.N + 1)
    alphas = Vector{Vector{Float64}}()
    activos = Vector{BitVector}()
    topes = Vector{BitVector}()
    σ = Float64[]
    el = Float64[]
    dis = Float64[]
    its = Int[]
    kkt = 0.0
    fijos = nodos_fijos(p)
    for d in carga
        previo = A
        r = alternancia_at2(previo, p, d; tol, iter)
        A = r.alpha
        u = respuesta_u(A, p, d)
        push!(alphas, copy(A))
        push!(activos, BitVector(A[i] <= previo[i] && !fijos[i] && A[i] < p.amax
                                 for i in eachindex(A)))
        push!(topes, BitVector(A[i] >= p.amax for i in eachindex(A)))
        push!(its, r.iteraciones)
        kkt = max(kkt, r.kkt)
        push!(σ, u.σ)
        push!(el, u.elastica)
        push!(dis, energia_disipada(A, p))
    end
    (alpha = alphas, sigma = σ, elastica = el, disipada = dis, total = el .+ dis,
     carga = collect(carga), iteraciones = its, activos = activos, tope = topes, fijos = fijos,
     p = p, kkt = kkt,
     longitud = [longitud_danada(a, p) for a in alphas])
end

# Rama homogénea: la solución cerrada contra la que se compara el solver.

# ε_c = √(G_c/(Eℓ)), la escala de deformación del AT2.
eps_c(p::ParAT2) = sqrt(p.Gc / (p.E * p.ℓ))
alpha_homogenea(ε::Real, p::ParAT2) = ε^2 / (ε^2 + eps_c(p)^2)
sigma_homogenea(ε::Real, p::ParAT2) = p.E * ε * (1 - alpha_homogenea(ε, p))^2

# El pico de σ(ε) está en ε = ε_c/√3 y vale (3√3/16)√(E G_c/ℓ): la resistencia va como ℓ^{-1/2}.
const AT2_PICO = 3 * sqrt(3) / 16
sigma_pico_homogenea(p::ParAT2) = AT2_PICO * sqrt(p.E * p.Gc / p.ℓ)
eps_pico_homogenea(p::ParAT2) = eps_c(p) / sqrt(3)
