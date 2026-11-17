# F4 · Objetivo 5 · El daño por gradiente AT2 en 1D: el modelo de EDP al que se traslada el control.
#
# El objetivo 5 de la propuesta dice que el control del caudal se extiende a un modelo de fractura por
# campo de fase. Esto es la cuenta que convierte «se trasladaría» en «acá está», con el MISMO modelo
# que la tesis usa para el daño del macizo (AT2, w(α) = α², sin umbral), de forma que el trabajo sirve
# en los dos lados. La energía es
#
#     E(u, α) = ∫ ½ (1-α)² E u'² dx  +  (G_c/2) ∫ ( α²/ℓ + ℓ α'² ) dx,
#
# que es AT2 con la normalización de Pham–Marigo: c_w = 4∫₀¹√(w) = 2 para w(α) = α². Dos cosas la hacen
# cómoda y las dos se usan abajo:
#
#  · En 1D sin carga de volumen la tensión es CONSTANTE en la barra, así que el paso en u sale en forma
#    cerrada —σ = E d / (h Σ_e (1-α_e)^{-2})— y no hay ningún sistema lineal que resolver. La energía
#    elástica reducida queda Ê_el(α; d) = E d² / (2h Σ_e (1-α_e)^{-2}).
#  · AT2 es CUADRÁTICA en α, así que el paso en α es un problema cuadrático convexo con cotas, sin
#    linealizar nada. La irreversibilidad α̇ ≥ 0 entra como cota inferior nodo por nodo (α ≥ α del paso
#    de carga anterior) y se impone por proyección dentro del propio Gauss–Seidel.
#
# El control del esquema es la rama HOMOGÉNEA, que es el único caso con solución cerrada: con cotas
# naturales en α, α uniforme es punto estacionario EXACTO también del problema discreto (se verifica
# nodo a nodo, borde incluido), así que el solver tiene que reproducir
#
#     α(ε) = ε²/(ε² + ε_c²),   σ(ε) = E ε (1-α)²,   ε_c² = G_c/(E ℓ),
#
# y su pico σ_c = (3√3/16)·√(E G_c/ℓ) en ε = ε_c/√3, a tolerancia de solver y no de discretización.
#
# Para que haya LOCALIZACIÓN —y entonces una carga de ruptura estructural que converger en malla— hay
# que romper la simetría: `defecto` baja G_c con una gaussiana de ancho FIJO (independiente de h), no
# con un elemento debilitado, que sería un defecto que se encoge al refinar.

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

# Sobre `fija_borde`, que NO es un detalle de implementación. Con cotas naturales en α el nodo del
# borde tiene la mitad de la masa disipativa de un nodo interior y un solo vecino que lo penalice en
# α', así que es el punto más débil de la barra: la localización se va al EXTREMO y no al defecto, por
# grande que sea el defecto. Medido el 14-nov con defecto = 0,3: α llegó a 1 en el nodo 1 mientras el
# centro quedaba en 0,27. Fijar α = 0 en los extremos es la condición estándar (extremos reforzados) y
# es la que usan las pruebas de localización. La rama homogénea, en cambio, EXIGE las cotas naturales:
# con α fijo en el borde no hay solución uniforme. Por eso el campo existe y no hay un único modo.

ParAT2(; L = 1.0, E = 1.0, Gc = 1.0, ℓ = 0.5, N = 60, defecto = 0.0, ancho = 0.15,
        amax = 1 - 1e-6, fija_borde = false) =
    ParAT2(L, E, Gc, ℓ, N, defecto, ancho, amax, fija_borde)

"Los nodos pinchados en α = 0: ninguno, o los dos extremos si `fija_borde`."
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

"G_c por elemento. Con `defecto = 0` es constante y la barra admite la solución homogénea exacta."
gc_elem(p::ParAT2) = p.defecto == 0 ? fill(p.Gc, p.N) :
    [p.Gc * (1 - p.defecto * exp(-((x - p.L / 2) / p.ancho)^2)) for x in centros_at2(p)]

"α en el punto medio de cada elemento, que es donde se evalúa la degradación elástica."
alpha_elem(A::AbstractVector) = (A[1:(end - 1)] .+ A[2:end]) ./ 2

"∫α dx por trapecios: la medida de «cuánto largo de barra está dañado» que usa el costo del control."
longitud_danada(A::AbstractVector, p::ParAT2) = paso_malla(p) * sum(alpha_elem(A))

# ---- el paso en u, en forma cerrada --------------------------------------------------------------

"""
    respuesta_u(A, p, d)

Minimización en u a α fijo, con u(0) = 0 y u(L) = d. Devuelve σ (constante en la barra), los saltos
de u por elemento, y S = Σ_e (1-α_e)^{-2} con las potencias de (1-α_e) que después usan el gradiente
y el Hessiano reducidos.
"""
function respuesta_u(A::AbstractVector, p::ParAT2, d::Real)
    αe = alpha_elem(A)
    w = 1 .- αe
    s = @. 1 / w^2
    S = sum(s)
    h = paso_malla(p)
    σ = p.E * d / (h * S)
    (σ = σ, S = S, s = s, αe = αe, w = w, du = d .* s ./ S, elastica = 0.5 * σ * d)
end

# ---- la parte disipada: energía, y su forma cuadrática ------------------------------------------

"""
    energia_disipada(A, p)

(G_c/2)∫(α²/ℓ + ℓα'²), con ∫α² exacta sobre P1 (matriz de masa, no punto medio) y α' constante por
elemento. `hessiano_disipado` es su Hessiano; que los dos coincidan es una prueba.
"""
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

"Hessiano de `energia_disipada`: masa (h/3ℓ, h/6ℓ) más rigidez (ℓ/h, -ℓ/h), ensamblado por elemento."
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

# ---- el paso en α: un problema cuadrático con cotas ---------------------------------------------

"""
    qp_alpha(u, p)

El paso en α a u fijo: min ½AᵗHA − bᵗA sobre α_prev ≤ A ≤ amax. H = H_disipado + H_elástico(u) y b
viene del término lineal de (1-α)²; `u` es la salida de `respuesta_u`.
"""
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

"""
    pgs!(A, H, b, lo, hi; iter, tol)

Gauss–Seidel proyectado para ½AᵗHA − bᵗA con cotas. Con H definida positiva cada barrido baja la
energía y converge; es exactamente minimizar coordenada por coordenada, así que la proyección es
`clamp` y no hace falta nada más.
"""
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

"""
    qp_cotas!(A, H, b, lo, hi; iter)

El mismo problema, por CONJUNTO ACTIVO primal: cada vuelta resuelve el sistema tridiagonal
restringido a los nodos libres —O(n), contra los O(n²) barridos que necesita Gauss–Seidel en una
malla fina—, manda a la cota los que la violan y libera los que tienen el multiplicador del signo
equivocado. Cierra con un pulido de `pgs!` y devuelve el residuo KKT: la prueba se hace sobre ESE
número, no sobre el método, así que si el conjunto activo ciclara se vería.
"""
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
        # Umbral del multiplicador: en una etapa que repite la carga anterior los nodos quedan en la
        # cota con multiplicador CERO, y liberar por un -1e-16 hace ciclar el conjunto activo sin
        # cambiar la solución. El umbral es relativo a la escala del lado derecho.
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
    # Si el conjunto activo AGOTÓ sus `iter` vueltas sin converger, salió ciclando y `A` no cumple
    # KKT: el pulido de 50 barridos no alcanza a arreglarlo y el residuo se publicaba como ~1,6e-2.
    # En ese caso —y solo en ese— el pulido pasa a ser el solver, con tope ancho. Medido el 15-nov:
    # de 75 921 QP, 66 786 agotan el conjunto activo, 211 se mueven más de 1e-6, y las cinco cifras
    # de `kkt_por_etapas` quedan todas en 1,4–1,7e-14 (antes: tres en 1,6e-2).
    pgs!(A, H, b, lo, hi; iter = vueltas == iter ? 200_000 : 50, tol = 1e-15)
    g = H * A .- b
    kkt = 0.0
    for i in 1:n
        # Un nodo con lo == hi está PINCHADO: no tiene condición de optimalidad, y medirlo como si la
        # tuviera daba un residuo de 0,4 con la solución exacta (medido el 14-nov con α = 0 en el
        # borde). El residuo tiene que hablar solo de los nodos que el problema puede mover.
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

"""
    alternancia_at2(A_prev, p, d; tol, iter)

Minimización alternada en u y α para UNA etapa de carga, con irreversibilidad por proyección sobre
`A_prev`. Devuelve el α convergido, los barridos gastados y la sucesión de energías reducidas
Ê(A^m) = min_u E(u, A^m), que tiene que ser no creciente: ése es el control del esquema.
"""
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

"""
    resolver_at2(carga, p; tol, iter)

Forward completo: una etapa por entrada de `carga` (desplazamiento impuesto en x = L), con la
irreversibilidad arrastrada entre etapas. Devuelve la historia de α, σ, energías, ∫α dx y, para el
adjunto, qué nodos quedaron CONGELADOS en cada etapa (`activos` = pegados a la cota de
irreversibilidad, `tope` = pegados a amax, que debería estar siempre vacío).
"""
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

# ---- la rama homogénea, que es el control analítico ---------------------------------------------

"ε_c = √(G_c/(Eℓ)), la única escala de deformación del AT2."
eps_c(p::ParAT2) = sqrt(p.Gc / (p.E * p.ℓ))
alpha_homogenea(ε::Real, p::ParAT2) = ε^2 / (ε^2 + eps_c(p)^2)
sigma_homogenea(ε::Real, p::ParAT2) = p.E * ε * (1 - alpha_homogenea(ε, p))^2

# El pico de σ(ε) está en ε = ε_c/√3 y vale (3√3/16)√(E G_c/ℓ). Es la «relación conocida para AT2»:
# la resistencia va como ℓ^{-1/2}, que es la ley que la figura y la prueba de invariancia miden.
const AT2_PICO = 3 * sqrt(3) / 16
sigma_pico_homogenea(p::ParAT2) = AT2_PICO * sqrt(p.E * p.Gc / p.ℓ)
eps_pico_homogenea(p::ParAT2) = eps_c(p) / sqrt(3)
