# La atribución del 1e-11 del cruce de Lyapunov sobre la nominal del óptimo (§4.4 del informe), medida
# dentro de la cadena para que sus cifras salgan como macros (\cfcru...) y no tipeadas a mano.
#
# Es la medición de docs/marco/ataque-s5/ (CRUCE-LYAPUNOV-MEDIDO.md; scripts medir/piso.jl, cruce_nudos.jl,
# escala.jl, estado_rk4.jl) reunida en una sola función: el mismo cruce que `cruce_lyapunov_nominal`, con el
# integrador escrito a mano para fijar reltol/abstol y los tstops desde afuera, y con el número de pasos de Ψ.
# Control de reproducción: el caso «óptimo, tal cual, reltol 1e-12» tiene que dar BIT A BIT el
# `error_rel_maximo` publicado en objetivo2_lineal.json; si no, se detiene (la medición no sería la del json).
#
# Uso suelto, desde la raíz:  julia/bin/jl --project=julia julia/cruce_atribucion.jl
# Dentro de la cadena: julia/regenerar_final.jl lo incluye y escribe final/resultados/cruce_lyapunov_atribucion.json.
import Fractura
using LinearAlgebra, StaticArrays, QuadGK, OrdinaryDiffEqVerner, JSON
# `F` global y const, como en medir/piso.jl: con `F = Fractura` local dentro de las funciones el cruce
# da otro número (2,79e-11 contra 2,68e-11, medido el 19-nov) y el control de reproducción falla.
const F = Fractura

function _integ_cruce(g, y0, tspan; rtol, atol, nudos = Float64[])
    kw = isempty(nudos) ? NamedTuple() : (tstops = nudos,)
    s = solve(ODEProblem{false}((y, _, t) -> g(y, t), y0, tspan), Vern9();
              reltol = rtol, abstol = atol, kw...)
    s.retcode == ReturnCode.Success || error("cruce_atribucion: integrador: $(s.retcode)")
    s
end

function _cruce(nom; T = nom.T, rtol, atol, nudos = Float64[])
    fW = function (Y, t)
        AB = F.matrices_lineales(nom, t); Wt = SMatrix{2, 2}(Y)
        vec(AB[1] * Wt + Wt * AB[1]' + AB[2] * AB[2]')
    end
    WL = F.simetrizar(SMatrix{2, 2}(_integ_cruce(fW, vec(0.0 * F.I2), (0.0, T); rtol, atol, nudos)(T)))
    fV = function (Y, t)
        AB = F.matrices_lineales(nom, t); Vt = SMatrix{2, 2}(Y)
        vec(-(AB[1]' * Vt + Vt * AB[1] + F.C_PRESION' * F.C_PRESION))
    end
    ML = F.simetrizar(SMatrix{2, 2}(_integ_cruce(fV, vec(0.0 * F.I2), (T, 0.0); rtol, atol, nudos)(0.0)))
    Ψs = _integ_cruce((Y, t) -> vec(-SMatrix{2, 2}(Y) * F.matrices_lineales(nom, t)[1]),
                      vec(F.I2), (T, 0.0); rtol, atol, nudos)
    Φs = _integ_cruce((Y, t) -> vec(F.matrices_lineales(nom, t)[1] * SMatrix{2, 2}(Y)),
                      vec(F.I2), (0.0, T); rtol, atol, nudos)
    Wq = SMatrix{2, 2}(quadgk(s -> (g = SMatrix{2, 2}(Ψs(s)) * F.∂fu(nom.x(s), nom.p); g * g'),
                              0.0, T; rtol = 1e-12, atol = 1e-14)[1])
    Mq = SMatrix{2, 2}(quadgk(s -> (g = F.C_PRESION * SMatrix{2, 2}(Φs(s)); g' * g),
                              0.0, T; rtol = 1e-12, atol = 1e-14)[1])
    λWq = eigvals(Symmetric(Matrix(Wq))); λMq = eigvals(Symmetric(Matrix(Mq)))
    eW = maximum(abs.(eigvals(Symmetric(Matrix(WL))) - λWq)) / max(1.0, maximum(abs, λWq))
    eM = maximum(abs.(eigvals(Symmetric(Matrix(ML))) - λMq)) / max(1.0, maximum(abs, λMq))
    (emax = max(eW, eM), pasos = length(Ψs.t), lambda_max_W = maximum(abs, λWq))
end

function cruce_atribucion()
    p = F.Par(); N = F.N_OBJ2
    q = F.optimizar(p, N).x
    nom = F.nominal_optima(p, q)                               # estado RK4-Hermite, la del json
    t = collect(F.malla(p, N)); nudos = t[2:(end - 1)]
    idx(τ) = clamp(searchsortedlast(t, τ), 1, length(q))
    nomV = F.nominal_control(p, τ -> q[idx(τ)], nom.x(0.0); T = p.T)   # misma ley, estado Vern9
    suave = F.nominales_semilla(p)["trayectoria_presurizada"]
    nudos_suave = collect(range(0.0, suave.T; length = N + 1))[2:(end - 1)]
    t12 = (rtol = 1e-12, atol = 1e-14); t14 = (rtol = 1e-14, atol = 1e-16)

    A   = _cruce(nom; t12...)
    B   = _cruce(nom; t12..., nudos)
    D   = _cruce(suave; t12...)
    E   = _cruce(suave; t12..., nudos = nudos_suave)
    A14 = _cruce(nom; t14...)
    V   = _cruce(nomV; t12...)

    publicado = JSON.parsefile(joinpath(F.RES_FINAL, "objetivo2_lineal.json"))["cruce_lyapunov_en_la_nominal"]["error_rel_maximo"]
    A.emax == publicado || error("cruce_atribucion: el caso A no reproduce el json publicado: $(A.emax) ≠ $publicado")

    Dict(
        "metodo" => "Cruce de Lyapunov contra cuadratura del propagador (el de cruce_lyapunov_en_la_nominal), " *
                    "con Vern9 a mano para fijar tolerancias y tstops. optimo = nominal del óptimo (RK4-Hermite, " *
                    "N_OBJ2 tramos); suave = trayectoria_presurizada de la semilla; tstops en los nudos interiores " *
                    "(en la suave, nudos arbitrarios equiespaciados, control negativo). El caso optimo/reltol 1e-12 " *
                    "reproduce bit a bit el error_rel_maximo de objetivo2_lineal.json (exigido).",
        "n_nudos_interiores" => length(nudos),
        "e_max" => Dict("optimo" => A.emax, "optimo_tstops" => B.emax, "suave" => D.emax,
                        "suave_tstops" => E.emax, "optimo_reltol_menos_catorce" => A14.emax,
                        "misma_ley_estado_vern" => V.emax),
        # Pasos del propagador Ψ (no de las EDO de Lyapunov). Con tstops y a reltol 1e-14 se guardan desde el 19-nov
        # (juez B1/M1): con nudos el óptimo BAJA de pasos, así que su ganancia no es refinamiento.
        "pasos" => Dict("optimo" => A.pasos, "suave" => D.pasos, "optimo_tstops" => B.pasos,
                        "suave_tstops" => E.pasos, "optimo_reltol_menos_catorce" => A14.pasos),
        "razon" => Dict("brecha_optimo_suave" => A.emax / D.emax,
                        "ganancia_tstops_optimo" => A.emax / B.emax,
                        "ganancia_tstops_suave" => D.emax / E.emax,
                        "ganancia_reltol" => A.emax / A14.emax,
                        "factor_estado" => A.emax / V.emax),
        "escala" => Dict("lambda_max_W_optimo" => A.lambda_max_W, "lambda_max_W_suave" => D.lambda_max_W,
                         "diferencia_pct" => 100 * abs(A.lambda_max_W / D.lambda_max_W - 1)),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    o = cruce_atribucion()
    destino = joinpath(Fractura.RES_FINAL, "cruce_lyapunov_atribucion.json")
    open(destino, "w") do io
        write(io, JSON.json(o; pretty = 2), "\n")
    end
    println("escrito: ", destino)
    println(JSON.json(o; pretty = 2))
end
