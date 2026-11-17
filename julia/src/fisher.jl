# F3.2 · La matriz de información de Fisher: cuán bien se separan los parámetros midiendo solo la
# presión del pozo, una vez que el calibre ya explicó qué NO se puede separar nunca.
#
#     M(q) = Σ_k s_k s_kᵀ / σ²,     s_k = ∂P(τ_k)/∂θ,
#
# con el modelo de observación del objetivo 2: en el pozo se mide SOLO la presión (C = [0 1]), con
# ruido gaussiano independiente de desviación σ. Cramér–Rao da cov(θ̂) ⪰ M⁻¹ para todo estimador no
# sesgado. Acá la cota se usa como PREDICCIÓN y se contrasta con la dispersión empírica de mínimos
# cuadrados sobre réplicas sintéticas: sin ese contraste, «la cota dice que se identifica» es una
# cuenta sobre el modelo, no un resultado sobre el estimador.
#
# Tres trampas, declaradas porque las tres dan números plausibles y falsos:
#
#  · Una varianza marginal chica CON LOS OTROS PARÁMETROS FIJOS no es identificabilidad. Lo que hay
#    que leer es la diagonal de M⁻¹, no 1/M_jj. El cociente inflacion_jj = sqrt(M⁻¹_jj · M_jj) ≥ 1
#    es exactamente cuánto se apoya θ_j en los demás; en el modelo linealizado y con ruido gaussiano
#    ese cociente ES el perfilado, en forma cerrada. El perfil no lineal lo aproxima el experimento
#    de réplicas, que no linealiza nada.
#  · Si la dispersión empírica queda muy por encima del ancho de la cota, la cota no mide nada en ese
#    parámetro: se reporta «no identificable con este diseño», y no «la cota falla».
#  · σ fija la escala de M pero NO los cocientes: la razón empírica/cota y la inflación son
#    invariantes en σ mientras el estimador siga en el régimen lineal. Por eso el experimento corre
#    con σ chico; con σ grande lo que se mide es la no linealidad, que es otra pregunta.

const SIGMA_FISHER = 1e-3

_media(v) = sum(v) / length(v)
_desv(v) = (m = _media(v); sqrt(sum(abs2, v .- m) / (length(v) - 1)))

"""
    fisher(p, ufun, ts; pars, σ, relativa)

M(q) para el diseño `ufun` medido en `ts`. Con `relativa = true` las derivadas son respecto de
log θ_j, así que M queda en unidades comparables entre parámetros de magnitudes distintas.
"""
fisher(p::Par, ufun, ts; pars = PARS_SENS, σ = SIGMA_FISHER, relativa = false) =
    (J = sens_presion(p, ufun, ts; pars, relativa); Symmetric(Matrix(J' * J) ./ σ^2))

cota_cramer_rao(M) = sqrt.(diag(inv(M)))

# sqrt(M⁻¹_jj · M_jj): 1 si θ_j no se confunde con ningún otro, grande si se apoya en ellos.
inflacion(M) = sqrt.(diag(inv(M)) .* diag(M))

"""
    ajustar(p, ufun, ts, datos; pars)

Mínimos cuadrados sobre `pars` por Gauss–Newton, usando las MISMAS sensibilidades de la EDO
aumentada: θ ← θ + (JᵗJ)⁻¹Jᵗr. Devuelve el vector ajustado y cuántas iteraciones gastó.
"""
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

"""
    experimento_cramer_rao(p, ufun, ts; pars, σ, replicas, semilla)

Genera `replicas` conjuntos de datos sintéticos con el θ verdadero y ruido σ, los ajusta por
mínimos cuadrados, y compara la desviación empírica de los estimadores con sqrt(diag(M⁻¹)).
"""
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
