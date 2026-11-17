# F3 · Sensibilidades paramétricas, y el calibre que explica por qué el problema inverso es singular.
#
# Dos piezas, que solo juntas contestan «¿se pueden identificar los parámetros midiendo la presión?».
#
#   1 · sensibilidades / sens_presion     ∂P(τ_k)/∂θ en el modelo ADIMENSIONAL, con θ entre
#       (n, c, η, P0), por EDO aumentada (variacional): el estado y sus derivadas respecto de los
#       parámetros se integran JUNTOS, con el mismo Vern9 y las mismas tolerancias del resto. Las
#       diferencias finitas centradas son el CONTROL de esto, no el método (prueba 3.1).
#
#   2 · presion_dim / calibre_* / invariantes     el modelo DIMENSIONAL de seis parámetros
#       (E', K_Ic, λ, v₀, C, L₀) —el del `.py` de MA5306— y el grupo de DOS parámetros que deja la
#       presión medida exactamente igual.
#
# Por qué la pieza 2 es el resultado y no un detalle: el `.py` de MA5306 encontró por SVD que el
# Jacobiano de la presión respecto de los seis parámetros tiene dos valores singulares ~1e-10
# relativos, y predijo el núcleo. Un valor singular chico, por sí solo, no distingue «dirección
# exactamente plana» de «dirección mal condicionada»: la SVD de un Jacobiano por diferencias finitas
# no puede bajar de su propio piso de ruido. Acá el grupo está DERIVADO del sistema y verificado a
# amplitud finita —a = 2, b = 3, no infinitesimal—, que es donde las dos cosas se separan.
#
# La derivación, que es tres líneas. Con L̃ = aL, p̃ = p y t̃ = t, el sistema dimensional
#
#     L̇ = v₀ (p√(πL)/K)ⁿ,      ṗ = (q − 2λLp − (4π/E) p L L̇)/(C + 2πL²/E)
#
# vuelve a tener la misma forma si E → a²E, λ → λ/a, L₀ → aL₀ y v₀ → a^{1−n/2} v₀, con C y q
# intactos. Y K y v₀ entran SOLO por v₀/Kⁿ, así que K → bK con v₀ → bⁿ v₀ tampoco se nota. Son dos
# direcciones, y la presión medida p(t_k) no cambia en ninguna: el problema inverso de seis
# parámetros es singular por construcción, con corrango exactamente 2.
#
# Lo identificable son los cuatro invariantes del grupo: C, E/L₀², λL₀ y v₀L₀^{n/2−1}/Kⁿ. Esa es la
# lectura para el póster, y no es pesimista: dice QUÉ medir (cuatro combinaciones, no seis números) y
# dice qué hay que fijar aparte —L₀ y K, o v₀— para que los seis se separen. El `.py` lo muestra:
# con L₀ y v₀ conocidos, el valor singular más chico sube de 2e-10 a 2e-5.
#
# Lo que esto NO dice: nada sobre la identificabilidad PRÁCTICA de las combinaciones que sí se
# identifican. Eso es la matriz de Fisher y la frontera, y ahí la pregunta es cuantitativa.

# ---- 1 · sensibilidades del modelo adimensional -------------------------------------------------

const PARS_SENS = (:n, :c, :η, :P0)

# Cuántas veces se evaluó ∂fθ fuera del dominio físico (ℓ, P > 0). En las ETAPAS INTERNAS de un paso
# tentativo de Vern9 el estado puede salirse aunque la trayectoria aceptada no lo haga: medido el
# 14-nov en el barrido de F3.3, con ℓ = 1,03 y P = −2,1e-3 en un paso que el controlador de error
# después rechazó. Abortar ahí mata la integración por un estado que no se usa; evaluar en el borde
# del dominio y CONTARLO deja la señal intacta. Las pruebas exigen 0 en el camino normal.
const PISO_SENS = Ref(0)
const EPS_SENS = 1e-12

# ∂f/∂θ_j. P0 y ℓ0 no aparecen en f: su sensibilidad viene toda de la condición inicial y acá dan 0.
function ∂fθ(x, u, p::Par, par::Symbol)
    ℓ, P, n = x[1], x[2], p.n
    d = D(x, p)
    if par === :n
        if !(ℓ > 0 && P > 0)
            PISO_SENS[] += 1
            ℓ = max(ℓ, EPS_SENS)
            P = max(P, EPS_SENS)
        end
        F = ℓ^(n / 2) * P^n
        Fn = F * (log(ℓ) / 2 + log(P))            # ∂F/∂n, con F = ℓ^{n/2}Pⁿ
        SVector(Fn, -2ℓ * P * Fn / d)
    elseif par === :c
        F = ℓ^(n / 2) * P^n
        SVector(0.0, -(u - p.η * ℓ * P - 2ℓ * P * F) / d^2)   # c entra solo por D = c + ℓ²
    elseif par === :η
        SVector(0.0, -ℓ * P / d)
    elseif par === :ℓ0 || par === :P0
        SVector(0.0, 0.0)
    else
        error("∂fθ: parámetro $par no implementado")
    end
end

# ∂x(0)/∂θ_j: cero salvo para las condiciones iniciales, que son el parámetro mismo.
_s0(par::Symbol) = par === :ℓ0 ? SVector(1.0, 0.0) : par === :P0 ? SVector(0.0, 1.0) : SVector(0.0, 0.0)

"""
    sensibilidades(p, ufun, ts; pars = PARS_SENS)

Para cada t de `ts`, el estado y la matriz S = ∂x(t)/∂θ (2 × length(pars)), integrando la EDO
aumentada Ṡ_j = A(t)S_j + ∂f/∂θ_j junto al estado. `ufun` es el control como función del tiempo.
"""
function sensibilidades(p::Par, ufun, ts; pars = PARS_SENS)
    m = length(pars)
    y0 = vcat(SVector(p.ℓ0, p.P0), reduce(vcat, map(_s0, pars)))
    s = integrar(y0, (0.0, maximum(ts))) do y, t
        x = SVector(y[1], y[2])
        u = ufun(t)
        A = ∂f(x, u, p)
        cols = ntuple(j -> A * SVector(y[1 + 2j], y[2 + 2j]) + ∂fθ(x, u, p, pars[j]), m)
        vcat(f(x, u, p), reduce(vcat, cols))
    end
    [(y = s(t); (x = SVector(y[1], y[2]), S = SMatrix{2, m}(y[3:end]))) for t in ts]
end

"""
    sens_presion(p, ufun, ts; pars = PARS_SENS, relativa = false)

∂P(t_k)/∂θ_j como matriz de length(ts) × length(pars): una fila por instante de medición, una
columna por parámetro. Con `relativa = true`, ∂P/∂log θ_j = θ_j ∂P/∂θ_j, que es lo comparable entre
parámetros de magnitudes distintas y lo que usa el `.py` de MA5306.
"""
function sens_presion(p::Par, ufun, ts; pars = PARS_SENS, relativa = false)
    r = sensibilidades(p, ufun, ts; pars)
    J = [r[k].S[2, j] for k in eachindex(r), j in eachindex(pars)]
    relativa ? J .* reshape([Float64(getfield(p, par)) for par in pars], 1, :) : J
end

# La presión del modelo adimensional, por integración directa: el control de las sensibilidades por
# diferencias finitas sale de acá, y no de la EDO aumentada.
function presion_adim(p::Par, ufun, ts)
    s = integrar((y, t) -> f(y, ufun(t), p), SVector(p.ℓ0, p.P0), (0.0, maximum(ts)))
    [s(t)[2] for t in ts]
end

# ---- 2 · el modelo dimensional de seis parámetros, y su calibre ---------------------------------

# Los mismos valores base del `.py` de MA5306 (y p0 = 0,3 dimensional, q ≡ 1, T = 4).
const DIM_BASE = (E = 1.0, K = 1.0, λ = 0.2, v0 = 1.0, C = 1.0, L0 = 1.0)

"""
    presion_dim(θ, n, p0, qfun, ts)

p(t_k) del modelo DIMENSIONAL con θ = (E, K, λ, v0, C, L0). Integrado con Vern9 y las tolerancias de
`integrar`; el `.py` de MA5306 usa Radau, así que comparar los dos es cruzar dos implementaciones.
"""
function presion_dim(θ, n, p0, qfun, ts)
    s = integrar(SVector(θ.L0, p0), (0.0, maximum(ts))) do y, t
        L, P = y[1], max(y[2], 0.0)
        Ld = θ.v0 * (P * sqrt(π * L) / θ.K)^n
        SVector(Ld, (qfun(t) - 2θ.λ * L * P - 4π / θ.E * P * L * Ld) / (θ.C + 2π * L^2 / θ.E))
    end
    [s(t)[2] for t in ts]
end

# Las dos direcciones del calibre, en perturbaciones LOGARÍTMICAS de (E, K, λ, v₀, C, L₀).
calibre_escala(n) = (E = 2.0, K = 0.0, λ = -1.0, v0 = 1 - n / 2, C = 0.0, L0 = 1.0)
calibre_vk(n) = (E = 0.0, K = 1.0, λ = 0.0, v0 = Float64(n), C = 0.0, L0 = 0.0)

# θ movido una cantidad `s` a lo largo de una dirección logarítmica: θ_k → θ_k exp(s·dir_k).
mover_log(θ, dir, s) = NamedTuple{keys(θ)}(map(k -> getfield(θ, k) * exp(s * getfield(dir, k)), keys(θ)))

# Los CUATRO invariantes del calibre: las combinaciones que la presión medida sí determina.
invariantes(n) = ("C" => (E = 0.0, K = 0.0, λ = 0.0, v0 = 0.0, C = 1.0, L0 = 0.0),
                  "E/L0^2" => (E = 1.0, K = 0.0, λ = 0.0, v0 = 0.0, C = 0.0, L0 = -2.0),
                  "lambda*L0" => (E = 0.0, K = 0.0, λ = 1.0, v0 = 0.0, C = 0.0, L0 = 1.0),
                  "v0*L0^(n/2-1)/K^n" => (E = 0.0, K = -Float64(n), λ = 0.0, v0 = 1.0, C = 0.0,
                                          L0 = n / 2 - 1))

monomio(θ, w) = prod(getfield(θ, k)^getfield(w, k) for k in keys(θ))

"""
    jacobiano_dim(θ, n, p0, qfun, ts; h = 1e-5, libres = keys(θ))

∂p(t_k)/∂log θ_j por diferencias finitas logarítmicas centradas, igual que el `.py` de MA5306: una
fila por instante, una columna por parámetro libre. Acá las diferencias finitas son el método, no el
control, porque lo que se mide es el núcleo y no el valor de la derivada.
"""
function jacobiano_dim(θ, n, p0, qfun, ts; h = 1e-5, libres = keys(θ))
    cols = map(libres) do k
        d = NamedTuple{(k,)}((1.0,))
        (presion_dim(mover_log(θ, merge(map(_ -> 0.0, θ), d), h), n, p0, qfun, ts) -
         presion_dim(mover_log(θ, merge(map(_ -> 0.0, θ), d), -h), n, p0, qfun, ts)) / 2h
    end
    hcat(cols...)
end

"""
    residuo_nucleo(J, dirs)

Cuánto se sale el núcleo numérico de J (sus `length(dirs)` últimas direcciones singulares derechas)
del subespacio que generan `dirs`. Es la métrica que imprime el `.py` de MA5306.
"""
function residuo_nucleo(J, dirs)
    P = hcat((collect(Float64, values(d)) for d in dirs)...)
    Q = Matrix(qr(P).Q)[:, 1:size(P, 2)]
    N = svd(J).V[:, end-length(dirs)+1:end]
    norm(N - Q * (Q' * N))
end
