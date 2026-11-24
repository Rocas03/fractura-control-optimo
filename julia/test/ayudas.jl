using Test, LinearAlgebra, QuadGK, Random

# Datos y funciones auxiliares de las pruebas.
@isdefined(SEMILLA) || (const SEMILLA = [0.8161, 0.4346])
const PAR = Par()
# Costos del método directo en el caso base, por N, como referencia de regresión.
const J_DIRECTO = Dict(25 => 0.37189675845351483, 50 => 0.3718706880094072,
                       100 => 0.3718641717379027, 200 => 0.37186254275375014)
# Dos vectores pseudoaleatorios fijos de 16 componentes, para controles y direcciones de prueba reproducibles.
const Z4703_A = [1.1101498731268962, -0.7430743313643473, -0.1889457562947505, 1.762304604308583,
                 0.06515654806700127, -1.4963277159880115, 1.4814237006854474, 0.9790828129816159,
                 1.3172588389435047, 0.47581029395806695, 0.5401358771577003, 0.555811537951228,
                 -0.033306911585706256, 1.1670923063244798, -0.055640185238523786, 0.5168164656944366]
const Z4703_B = [0.9593140929800286, -0.6010687870401767, -0.051765491716739144, 0.01819739677957788,
                 0.5944612495785925, 0.5539033658093443, 0.3196320224666674, 0.6027204799683691,
                 0.0831795068614417, 0.6236807011659488, 0.20464517094567267, 1.5152688327613608,
                 1.214494583443103, -0.9533762357061356, 1.5788847730487372, -0.24205926254636012]

allclose(a, b; rtol = 1e-5, atol = 1e-8) = all(abs.(a .- b) .<= atol .+ rtol .* abs.(b))

# El método directo refinando N, cada nivel desde el óptimo del anterior; devuelve el último.
function cascada(par, Ns)
    q = nothing
    res = nothing
    for N in Ns
        res = optimizar(par, N; q0 = q === nothing ? zeros(N) : repeat(q, inner = 2))
        q = res.x
    end
    res, q
end

# Problema LQ de doble integrador con solución exacta V = xᵀP(0)x / 2 por Riccati, para validar el
# esquema semi-Lagrangiano de HJB.
const T_LQ = 1.0
_lq_f(X, U) = (X[2], U)
_lq_costo(X, U) = 0.5 * (X[1]^2 + X[2]^2 + U^2)
_lq_g(X) = 0.5 * (X[1]^2 + X[2]^2)
_ejes(a, b, n) = (range(a, b, length = n), range(a, b, length = n))

function _riccati(t)
    A = [0.0 1.0; 0.0 0.0]
    B = reshape([0.0, 1.0], 2, 1)
    s = integrar((p, _) -> vec(-(A' * reshape(p, 2, 2) + reshape(p, 2, 2) * A -
                                  reshape(p, 2, 2) * B * B' * reshape(p, 2, 2) + I(2))),
                 vec(Matrix(1.0I, 2, 2)), (T_LQ, t))
    reshape(s.u[end], 2, 2)
end

# Error máximo de V(0) contra la solución de Riccati, en el interior de la malla para que no pese el borde.
function _error_lq(n, nt, nu; margen = 1.0)
    ejes = _ejes(-3.0, 3.0, n)
    V, _ = semi_lagrangiano(_lq_f, _lq_costo, _lq_g, ejes, range(-4.0, 4.0, length = nu), T_LQ, nt,
                            (:extrapolar, :extrapolar))
    P0 = _riccati(0.0)
    G0 = [x for x in ejes[1], y in ejes[2]]
    G1 = [y for x in ejes[1], y in ejes[2]]
    exacto = 0.5 .* (P0[1, 1] .* G0 .^ 2 .+ 2 .* P0[1, 2] .* G0 .* G1 .+ P0[2, 2] .* G1 .^ 2)
    interior = (abs.(G0) .<= margen) .& (abs.(G1) .<= margen)
    maximum(abs.(V[1] .- exacto)[interior])
end

# Oscilador armónico con costo cuadrático: otro LQ con Riccati conocido.
const A_OSC = [0.0 1.0; -1.0 0.0]
const B_OSC = [0.0, 1.0]
_osc_f(x, u) = A_OSC * x + B_OSC * u
_osc_costo(x, u) = 0.5 * (dot(x, x) + u^2)
_osc_g(x) = 0.5 * dot(x, x)

function _riccati_osc(T)
    integrar((p, _) -> (P = reshape(p, 2, 2); vec(-(A_OSC' * P + P * A_OSC - (P * B_OSC) * (B_OSC' * P) + I(2)))),
             vec(Matrix(1.0I, 2, 2)), (T, 0.0))
end
