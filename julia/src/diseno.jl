# F3.3 · La frontera control ↔ información: cuánto costo de control hay que pagar para que la misma
# inyección, además de fracturar bien, permita medir los parámetros del macizo.
#
#     min_q  Ψ_γ(q) = J(q) − γ · log det M(q),      q constante a trozos en [0, qmax],
#
# con J el costo del objetivo 3 (el mismo `avance`, el mismo Ipopt) y M la matriz de Fisher de F3.2
# en unidades RELATIVAS (derivadas respecto de log θ): sin eso, log det M mezcla parámetros de
# magnitudes distintas y el criterio D-óptimo premia al que tiene el valor más grande.
#
# γ = 0 es el control óptimo puro; γ → ∞ es el experimento óptimo sin importar el costo. La frontera
# es el conjunto de óptimos intermedios, y lo que se lee en ella es el precio de la información.
#
# El gradiente, y la desviación declarada: ∇J es el ADJUNTO exacto que ya existe (`gradiente`); el del
# término de información es por DIFERENCIAS FINITAS centradas en q. El gradiente exacto pediría el
# adjunto del sistema aumentado (estado + sensibilidades), o sea las segundas derivadas ∂²f/∂x∂θ, y eso
# no cabe antes del corte del lunes. El costo de la decisión está acotado y medido: 2N integraciones
# aumentadas por gradiente (N = 8 tramos → ~0,1 s), tolerancia de Ipopt relajada a 1e-8 (con 1e-12 la
# convergencia la fija el error de las diferencias finitas, no el problema), y la prueba de consistencia
# direccional del gradiente ensamblado, que es la que caza un error de índice, de signo o de γ.

const NS_DISENO = 8
const TS_DISENO = collect(range(0.2, 2.0, length = 40))
const IPOPT_DISENO = (print_level = 0, tol = 1e-8, limited_memory_max_history = 30, max_iter = 500,
                      honor_original_bounds = "yes")

"""
    escalon(q, p)

El control constante a trozos de `q` sobre `malla(p, length(q))`, como función del tiempo: es el mismo
control que usa el método directo, pero evaluable en cualquier instante (el directo lo usa tramo a
tramo). La prueba 3.3 de consistencia de integradores compara las dos vías.
"""
function escalon(q, p)
    t = malla(p, length(q))
    τ -> q[clamp(searchsortedlast(t, τ), 1, length(q))]
end

# Ipopt evalúa el objetivo en puntos que pueden salir de la caja [0, qmax] (fase de restauración y
# relajación de cotas), y ahí el modelo de sensibilidades no está definido: con q < 0 la presión se va
# por debajo de cero y ∂fθ, que pide log(P√(πL)/K), aborta. El término de información se evalúa
# entonces en la PROYECCIÓN de q sobre la caja, que dentro de ella no cambia nada y fuera lo deja
# constante (gradiente nulo en esa dirección). `Q_PROYECTADO` cuenta cuántas veces pasó: si fuera en el
# óptimo, el resultado no valdría, y por eso se reporta.
const Q_PROYECTADO = Ref(0)

const Q_ULTIMO = Ref{Vector{Float64}}(Float64[])

function caja(q, p)
    any(x -> x < 0 || x > p.qmax, q) && (Q_PROYECTADO[] += 1)
    Q_ULTIMO[] = collect(float.(q))
    clamp.(q, 0.0, p.qmax)
end

fisher_diseno(q, p; ts = TS_DISENO, σ = SIGMA_FISHER, pars = PARS_SENS) =
    fisher(p, escalon(caja(q, p), p), ts; pars, σ, relativa = true)

"""
    logdet_fisher(q, p; ts, σ, pars)

log det M(q). Por Cholesky, que además certifica que M es definida positiva: si no lo es, el diseño no
identifica y el logaritmo no tiene sentido, así que se aborta en vez de devolver un número.
"""
function logdet_fisher(q, p; ts = TS_DISENO, σ = SIGMA_FISHER, pars = PARS_SENS)
    F = cholesky(Symmetric(Matrix(fisher_diseno(q, p; ts, σ, pars))), check = false)
    issuccess(F) || error("M(q) no es definida positiva en q = $q: el diseño no identifica")
    logdet(F)
end

"""
    grad_logdet_fisher(q, p; h, ...)

∇_q log det M(q) por diferencias finitas centradas, con paso relativo `h` y diferencia hacia adelante
en el borde q_k = 0 (un q negativo no es un caudal y el modelo no lo admite).
"""
function grad_logdet_fisher(q, p; ts = TS_DISENO, σ = SIGMA_FISHER, pars = PARS_SENS, h = 1e-5)
    g = similar(q)
    for k in eachindex(q)
        hk = h * max(1.0, abs(q[k]))
        mas = copy(q); mas[k] += hk
        if q[k] - hk < 0
            g[k] = (logdet_fisher(mas, p; ts, σ, pars) - logdet_fisher(q, p; ts, σ, pars)) / hk
        else
            menos = copy(q); menos[k] -= hk
            g[k] = (logdet_fisher(mas, p; ts, σ, pars) -
                    logdet_fisher(menos, p; ts, σ, pars)) / 2hk
        end
    end
    g
end

psi(q, p, γ; kw...) = avance(q, p)[1] - (γ == 0 ? 0.0 : γ * logdet_fisher(q, p; kw...))
grad_psi(q, p, γ; kw...) =
    gradiente(q, p) - (γ == 0 ? zero(q) : γ .* grad_logdet_fisher(q, p; kw...))

struct Diseno <: AbstractNLPModel{Float64, Vector{Float64}}
    meta::NLPModelMeta{Float64, Vector{Float64}}
    counters::Counters
    p::Par
    γ::Float64
    ts::Vector{Float64}
    σ::Float64
end

Diseno(p, γ, q0; ts = TS_DISENO, σ = SIGMA_FISHER) =
    Diseno(NLPModelMeta(length(q0); x0 = float.(q0), lvar = zeros(length(q0)),
                        uvar = fill(p.qmax, length(q0)), hess_available = false),
           Counters(), p, γ, collect(ts), σ)

NLPModels.obj(m::Diseno, q::AbstractVector) = psi(collect(q), m.p, m.γ; ts = m.ts, σ = m.σ)
NLPModels.grad!(m::Diseno, q::AbstractVector, g::AbstractVector) =
    (g .= grad_psi(collect(q), m.p, m.γ; ts = m.ts, σ = m.σ))

"""
    grad_proyectado(q, g, p)

El gradiente proyectado sobre el cono admisible: en los tramos interiores es g, y en los que están
pegados a una cota solo cuenta la componente que empujaría hacia adentro. Su norma es la medida de
optimalidad que vale acá — Ipopt reporta `:first_order` contra su propia tolerancia, y con un gradiente
por diferencias finitas esa tolerancia no se alcanza aunque el punto sea el óptimo.
"""
function grad_proyectado(q, g, p)
    [(q[k] ≤ 1e-9 && g[k] > 0) || (q[k] ≥ p.qmax - 1e-9 && g[k] < 0) ? 0.0 : g[k] for k in eachindex(q)]
end

"""
    optimizar_diseno(p, γ; q0, ts, σ)

Un punto de la frontera: el minimizador de Ψ_γ y lo que hay que reportar de él.
"""
function optimizar_diseno(p, γ; q0 = fill(0.5, NS_DISENO), ts = TS_DISENO, σ = SIGMA_FISHER, kw...)
    st = ipopt(Diseno(p, γ, q0; ts, σ); IPOPT_DISENO..., kw...)
    q = clamp.(st.solution, 0.0, p.qmax)
    # El piso del dominio de ∂fθ se activa en etapas internas de pasos rechazados durante la
    # optimización. Lo que no puede pasar es que se active en el PUNTO devuelto: eso se mide acá.
    PISO_SENS[] = 0
    logdet_fisher(q, p; ts, σ)
    piso_en_el_optimo = PISO_SENS[]
    M = Matrix(fisher_diseno(q, p; ts, σ))
    (; γ, q, J = avance(q, p)[1], logdetM = logdet_fisher(q, p; ts, σ),
       lambda_min = minimum(eigvals(Symmetric(M))), cota = cota_cramer_rao(Symmetric(M)),
       psi = psi(q, p, γ; ts, σ), exito = st.status == :first_order, status = st.status, iter = st.iter,
       piso_en_el_optimo, grad_proy = norm(grad_proyectado(q, grad_psi(q, p, γ; ts, σ), p)))
end

"""
    frontera(γs; p, ts, σ)

La frontera completa, por continuación: cada γ arranca del óptimo del γ anterior, que es lo que evita
que Ipopt caiga en mínimos distintos y rompa la monotonía por una razón numérica y no física.
"""
function frontera(γs; p = Par(), ts = TS_DISENO, σ = SIGMA_FISHER, q0 = fill(0.5, NS_DISENO))
    out = []
    q = collect(float.(q0))
    for γ in γs
        r = optimizar_diseno(p, γ; q0 = q, ts, σ)
        push!(out, r)
        r.exito && (q = copy(r.q))
    end
    identity.(out)
end

const GAMMAS_FRONTERA = (0.0, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0)

# El «codo» de la frontera: el γ en que el informe señala el tramo barato de la curva. Es una ELECCIÓN
# de lectura, no un resultado del cálculo, y por eso va declarada una sola vez acá y la usan tanto el
# json (bloque `codo` de `diseno_fisher`) como la figura (`GAMMAS_ROTULADAS` en figuras.jl).
const GAMMA_CODO = 3e-2

"""
    diseno_fisher(p; γs, ts, σ)

La frontera completa como resultado publicable: un punto por γ, sus criterios, y la LECTURA del
resultado (qué tan ciego es el óptimo de control), separada de los criterios porque la lectura puede
salir en cualquiera de los dos sentidos y el criterio no.
"""
function diseno_fisher(p = Par(); γs = GAMMAS_FRONTERA, ts = TS_DISENO, σ = SIGMA_FISHER)
    ref = optimizar(p, NS_DISENO; q0 = fill(0.5, NS_DISENO))
    fr = frontera(γs; p, ts, σ)
    θ = [Float64(getfield(p, k)) for k in PARS_SENS]
    pts = [Dict("gamma" => r.γ, "J" => r.J, "logdetM" => r.logdetM, "lambda_min" => r.lambda_min,
                "cota_relativa" => r.cota, "q" => r.q, "status" => string(r.status),
                "iter" => r.iter, "grad_proy" => r.grad_proy, "piso" => r.piso_en_el_optimo) for r in fr]
    J = [r.J for r in fr]
    ld = [r.logdetM for r in fr]
    lm = [r.lambda_min for r in fr]
    # El índice del parámetro peor identificado en γ = 0, que es el que la frontera tiene que arreglar.
    # `cota` YA ES RELATIVA: M se construye con derivadas respecto de log θ (`relativa = true` en
    # `fisher_diseno`), así que sqrt(diag(M⁻¹))_j es la desviación estándar de log θ_j, o sea el error
    # relativo de θ_j. Hasta el 15-nov-2026 esto se dividía OTRA VEZ por θ_j: eso escalaba cada
    # parámetro por su propia magnitud y rompía justo la comparación entre parámetros que el diseño en
    # log θ busca. No se divide. Medido al sacar la división: `peor` no cambia (sigue siendo η).
    peor = argmax(fr[1].cota)
    criterios = Dict(
        "gamma0_es_el_optimo_de_control" => abs(fr[1].J - ref.fun) < 1e-6 * max(1.0, abs(ref.fun)),
        "J_no_decrece_con_gamma" => issorted(J),
        "logdetM_no_decrece_con_gamma" => issorted(ld),
        "lambda_min_no_decrece_con_gamma" => issorted(lm),
        "el_piso_del_dominio_no_se_activa_en_ningun_optimo" => all(r.piso_en_el_optimo == 0 for r in fr),
        "la_frontera_satura" => (ld[end] - ld[end-1]) < 0.05 * (ld[2] - ld[1]),
        "el_peor_parametro_mejora_al_menos_3x" =>
            (fr[1].cota[peor] / fr[end].cota[peor]) > 3.0)
    # Lectura, NO criterio: la frase «el óptimo de control es casi ciego» se escribe solo si λ_min en
    # γ = 0 es al menos dos órdenes de magnitud menor que en el extremo de diseño (plan, tarea 3.3).
    lectura = Dict("peor_parametro" => string(PARS_SENS[peor]),
                   "cota_relativa_peor_en_gamma0" => fr[1].cota[peor],
                   "cota_relativa_peor_en_gamma_max" => fr[end].cota[peor],
                   "factor_lambda_min" => lm[end] / lm[1],
                   "precio_en_J_del_extremo" => fr[end].J / fr[1].J - 1,
                   "optimo_de_control_casi_ciego" => lm[1] < lm[end] / 100,
                   "J_directo_referencia" => ref.fun)
    # El codo, con sus cifras: es lo que el informe necesita para decir «el tramo barato es el inicial»
    # sin tipear ningún número. `precio_en_J` y `mejora_peor` son cocientes contra γ = 0.
    i_codo = findfirst(r -> r.γ == GAMMA_CODO, fr)
    isnothing(i_codo) && error("GAMMA_CODO = $GAMMA_CODO no está en γs = $γs")
    codo = Dict("gamma" => fr[i_codo].γ, "J" => fr[i_codo].J,
                "cota_relativa_peor" => fr[i_codo].cota[peor],
                "precio_en_J" => fr[i_codo].J / fr[1].J - 1,
                "mejora_peor" => fr[1].cota[peor] / fr[i_codo].cota[peor])
    Dict("frontera" => pts, "criterios" => criterios, "lectura" => lectura, "codo" => codo,
         "diseno" => Dict("N_tramos" => NS_DISENO, "instantes" => length(ts), "sigma" => σ,
                          "pars" => string.(collect(PARS_SENS)), "theta" => θ,
                          "gradiente_informacion" => "diferencias finitas centradas, h = 1e-5 relativo",
                          "ipopt" => Dict(string(k) => string(v) for (k, v) in pairs(IPOPT_DISENO))))
end
