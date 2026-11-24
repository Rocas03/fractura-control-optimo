# Matriz de información de Fisher M(q) = Σ_k s_k s_kᵀ / σ², con s_k = ∂P(τ_k)/∂θ: cuán bien se separan los
# parámetros midiendo solo la presión del pozo, con ruido gaussiano independiente de desviación σ.
# Por Cramér–Rao, cov(θ̂) ⪰ M⁻¹ para todo estimador insesgado. La cota se contrasta con la dispersión de
# mínimos cuadrados sobre réplicas sintéticas; sin eso sería una cuenta sobre el modelo y no sobre el
# estimador.
# Lo que se lee es la diagonal de M⁻¹, no 1/M_jj: una varianza chica con los demás parámetros fijos no es
# identificabilidad. La inflación sqrt(M⁻¹_jj · M_jj) ≥ 1 mide cuánto se apoya θ_j en los otros.
# Si la dispersión empírica queda muy por encima de la cota, el parámetro no se identifica con ese diseño.
# σ fija la escala de M pero no los cocientes, mientras el estimador siga en régimen lineal; por eso el
# experimento usa σ chico.

const SIGMA_FISHER = 1e-3

_media(v) = sum(v) / length(v)
_desv(v) = (m = _media(v); sqrt(sum(abs2, v .- m) / (length(v) - 1)))

# M(q) para el diseño `ufun` medido en `ts`. Con `relativa = true` deriva respecto de log θ_j, para que M
# sea comparable entre parámetros de magnitudes distintas.
fisher(p::Par, ufun, ts; pars = PARS_SENS, σ = SIGMA_FISHER, relativa = false) =
    (J = sens_presion(p, ufun, ts; pars, relativa); Symmetric(Matrix(J' * J) ./ σ^2))

cota_cramer_rao(M) = sqrt.(diag(inv(M)))

# sqrt(M⁻¹_jj · M_jj): 1 si θ_j no se confunde con ningún otro, grande si se apoya en ellos.
inflacion(M) = sqrt.(diag(inv(M)) .* diag(M))

# Mínimos cuadrados sobre `pars` por Gauss–Newton, θ ← θ + (JᵗJ)⁻¹Jᵗr, con las mismas sensibilidades de
# la ecuación variacional. Devuelve θ ajustado y las iteraciones usadas.
function ajustar(p::Par, ufun, ts, datos; pars = PARS_SENS, iter = 20, tol = 1e-11)
    θ = [Float64(getfield(p, k)) for k in pars]
    usadas = iter
    for i in 1:iter
        q = con(p; NamedTuple{pars}(Tuple(θ))...)
        J = sens_presion(q, ufun, ts; pars)
        δ = (J' * J) \ (J' * (datos .- presion_adim(q, ufun, ts)))
        θ .+= δ
        if maximum(abs, δ ./ max.(abs.(θ), 1e-8)) < tol
            usadas = i
            break
        end
    end
    (θ = θ, iteraciones = usadas)
end

# Genera `replicas` conjuntos de datos sintéticos con el θ verdadero y ruido σ, los ajusta por mínimos
# cuadrados y compara la desviación empírica de los estimadores con la cota sqrt(diag(M⁻¹)).
function experimento_cramer_rao(p::Par, ufun, ts; pars = PARS_SENS, σ = SIGMA_FISHER,
                                replicas = 200, semilla = 20261114)
    M = fisher(p, ufun, ts; pars, σ)
    cota = cota_cramer_rao(M)
    base = presion_adim(p, ufun, ts)
    rng = MersenneTwister(semilla)
    est = zeros(replicas, length(pars))
    iters = zeros(Int, replicas)
    for i in 1:replicas
        r = ajustar(p, ufun, ts, base .+ σ .* randn(rng, length(ts)); pars)
        est[i, :] = r.θ
        iters[i] = r.iteraciones
    end
    θ = [Float64(getfield(p, k)) for k in pars]
    emp = [_desv(est[:, j]) for j in eachindex(pars)]
    sesgo = [_media(est[:, j]) for j in eachindex(pars)] .- θ
    (theta = θ, M = M, cota = cota, empirica = emp, razon = emp ./ cota,
     sesgo_sobre_cota = sesgo ./ cota, inflacion = inflacion(M),
     iteraciones_max = maximum(iters), replicas = replicas, sigma = σ)
end
