# Frontera entre costo de control e información: cuánto más hay que pagar en J para que la misma
# inyección, además de llevar la fractura al largo objetivo, permita estimar los parámetros.
# Se resuelve min_q Ψ_γ(q) = J(q) − γ log det M(q), con q constante a trozos en [0, qmax], J el costo del
# método directo y M la matriz de Fisher en unidades relativas (derivadas respecto de log θ); sin eso,
# log det M mezclaría parámetros de magnitudes distintas.
# γ = 0 es el control óptimo puro y γ grande es el experimento óptimo sin mirar el costo; la frontera son
# los óptimos intermedios y muestra el precio de la información.
# ∇J es el adjunto exacto de `gradiente`; el gradiente de log det M va por diferencias finitas centradas,
# porque el exacto pediría el adjunto del sistema con sensibilidades (segundas derivadas ∂²f/∂x∂θ). Por
# eso la tolerancia de Ipopt se relaja a 1e-8: más abajo manda el error de las diferencias finitas.

const NS_DISENO = 8
const TS_DISENO = collect(range(0.2, 2.0, length = 40))
const IPOPT_DISENO = (print_level = 0, tol = 1e-8, limited_memory_max_history = 30, max_iter = 500,
                      honor_original_bounds = "yes")

# El control constante a trozos de `q` sobre `malla(p, length(q))`, como función del tiempo: el mismo del
# método directo, pero evaluable en cualquier instante.
function escalon(q, p)
    t = malla(p, length(q))
    τ -> q[clamp(searchsortedlast(t, τ), 1, length(q))]
end

# Ipopt puede evaluar fuera de la caja [0, qmax], y con q < 0 la presión baja de cero y ∂fθ deja de estar
# definido. Por eso el término de información se evalúa en la proyección de q sobre la caja, y
# `Q_PROYECTADO` cuenta cuántas veces pasó; si pasara en el óptimo, el resultado no valdría.
const Q_PROYECTADO = Ref(0)

const Q_ULTIMO = Ref{Vector{Float64}}(Float64[])

function caja(q, p)
    any(x -> x < 0 || x > p.qmax, q) && (Q_PROYECTADO[] += 1)
    Q_ULTIMO[] = collect(float.(q))
    clamp.(q, 0.0, p.qmax)
end

fisher_diseno(q, p; ts = TS_DISENO, σ = SIGMA_FISHER, pars = PARS_SENS) =
    fisher(p, escalon(caja(q, p), p), ts; pars, σ, relativa = true)

# log det M(q) por Cholesky, que además confirma que M es definida positiva; si no lo es, el diseño no
# identifica y se aborta en vez de devolver un número.
function logdet_fisher(q, p; ts = TS_DISENO, σ = SIGMA_FISHER, pars = PARS_SENS)
    F = cholesky(Symmetric(Matrix(fisher_diseno(q, p; ts, σ, pars))), check = false)
    issuccess(F) || error("M(q) no es definida positiva en q = $q: el diseño no identifica")
    logdet(F)
end

# ∇_q log det M(q) por diferencias finitas centradas con paso relativo `h`, y hacia adelante en el borde
# q_k = 0, porque un caudal negativo no tiene sentido en el modelo.
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

# Gradiente proyectado sobre el cono admisible: g en los tramos interiores y, en los pegados a una cota,
# solo la componente que empuja hacia adentro. Su norma es la medida de optimalidad que se reporta, porque
# con un gradiente por diferencias finitas Ipopt no llega a su propia tolerancia aunque esté en el óptimo.
function grad_proyectado(q, g, p)
    [(q[k] ≤ 1e-9 && g[k] > 0) || (q[k] ≥ p.qmax - 1e-9 && g[k] < 0) ? 0.0 : g[k] for k in eachindex(q)]
end

# Un punto de la frontera: el minimizador de Ψ_γ y lo que se reporta de él.
function optimizar_diseno(p, γ; q0 = fill(0.5, NS_DISENO), ts = TS_DISENO, σ = SIGMA_FISHER, kw...)
    st = ipopt(Diseno(p, γ, q0; ts, σ); IPOPT_DISENO..., kw...)
    q = clamp.(st.solution, 0.0, p.qmax)
    # El piso de ∂fθ puede activarse en pasos rechazados durante la optimización, pero no en el punto
    # devuelto; eso se mide acá.
    PISO_SENS[] = 0
    logdet_fisher(q, p; ts, σ)
    piso_en_el_optimo = PISO_SENS[]
    M = Matrix(fisher_diseno(q, p; ts, σ))
    (; γ, q, J = avance(q, p)[1], logdetM = logdet_fisher(q, p; ts, σ),
       lambda_min = minimum(eigvals(Symmetric(M))), cota = cota_cramer_rao(Symmetric(M)),
       psi = psi(q, p, γ; ts, σ), exito = st.status == :first_order, status = st.status, iter = st.iter,
       piso_en_el_optimo, grad_proy = norm(grad_proyectado(q, grad_psi(q, p, γ; ts, σ), p)))
end

# La frontera por continuación: cada γ parte del óptimo del γ anterior, para que Ipopt no salte a otro
# mínimo y rompa la monotonía por una razón numérica.
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

# El codo de la frontera, donde el informe marca el fin del tramo barato. Es una elección de lectura, no
# un resultado del cálculo; la usan el json (bloque `codo`) y la figura (`GAMMAS_ROTULADAS` en figuras.jl).
const GAMMA_CODO = 3e-2

# La frontera completa como resultado: un punto por γ, sus criterios y la lectura (qué tan ciego es el
# óptimo de control), separada de los criterios porque puede salir en cualquiera de los dos sentidos.
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
    # Parámetro peor identificado en γ = 0, el que la frontera tiene que mejorar. `cota` ya es relativa: con
    # derivadas respecto de log θ, sqrt(diag(M⁻¹))_j es la desviación de log θ_j, o sea el error relativo
    # de θ_j, y no hay que dividir otra vez por θ_j.
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
    # Lectura, no criterio: «el óptimo de control es casi ciego» se afirma solo si λ_min en γ = 0 es al
    # menos cien veces menor que en el extremo de diseño.
    lectura = Dict("peor_parametro" => string(PARS_SENS[peor]),
                   "cota_relativa_peor_en_gamma0" => fr[1].cota[peor],
                   "cota_relativa_peor_en_gamma_max" => fr[end].cota[peor],
                   "factor_lambda_min" => lm[end] / lm[1],
                   "precio_en_J_del_extremo" => fr[end].J / fr[1].J - 1,
                   "optimo_de_control_casi_ciego" => lm[1] < lm[end] / 100,
                   "J_directo_referencia" => ref.fun)
    # El codo con sus cifras, para que el informe no tipee ningún número. `precio_en_J` y `mejora_peor` son
    # cocientes contra γ = 0.
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
