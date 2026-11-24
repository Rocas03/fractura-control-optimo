# Objetivo 1: la dinámica sin control y con caudal constante, en un barrido de n, η y U.
const NS_DINAMICA = (1, 2, 3)
const ETAS = (0.0, 0.2, 1.0)
const US = (0.0, 0.5, 1.0)

# Para cada caso se mide monotonía de ℓ, positividad de P, la cota P ≤ (S0 + Ut)/c y el balance de volumen
# S = S0 + Ut - η∫ℓP. Con η = 0 y U = 1 se comparan las pendientes log-log tardías de ℓ y P con las teóricas,
# y con P0 = 0 y U = 0 se comprueba que (ℓ0, 0) es equilibrio.
function objetivo1()
    base = Par()
    barrido, tendencias, equilibrios = Dict{String, Any}[], Dict{String, Any}[], Dict{String, Any}[]
    cotas_ok = true
    for n in NS_DINAMICA
        for η in ETAS
            p = con(base; n = float(n), η, T = 10.0)
            S0 = (p.c + p.ℓ0^2) * p.P0
            for U in US
                r = simular(p, U, 10.0, 401)
                fila = Dict{String, Any}("n" => n, "eta" => η, "U" => U, "ell_T" => r.ℓ[end], "P_T" => r.P[end],
                                         "S_T" => r.S[end], "max_P" => maximum(r.P), "min_P" => minimum(r.P),
                                         "min_diff_ell" => minimum(diff(r.ℓ)),
                                         "residuo_balance" => maximum(abs, r.S .- S0 .- (U .* r.t .- η .* r.I)))
                push!(barrido, fila)
                cotas_ok = cotas_ok && fila["min_diff_ell"] ≥ -1e-12 && fila["min_P"] ≥ -1e-12 &&
                           all(r.P .≤ (S0 .+ U .* r.t) ./ p.c .+ 1e-12)
            end
        end
        for T in (10.0, 1e4)
            r = simular(con(base; n = float(n), η = 0.0, T), 1.0, T, 11)
            pℓ = log(r.ℓ[end] / r.ℓ[2]) / log(10.0)
            pP = log(r.P[end] / r.P[2]) / log(10.0)
            tℓ = 2 * (n + 1) / (3n + 2)
            tP = -(n + 2) / (3n + 2)
            push!(tendencias, Dict{String, Any}("n" => n, "T" => T, "pendiente_ell" => pℓ, "teorica_ell" => tℓ,
                                                "pendiente_P" => pP, "teorica_P" => tP, "diferencia_ell" => pℓ - tℓ,
                                                "diferencia_P" => pP - tP))
        end
        for η in ETAS
            p = con(base; n = float(n), η, P0 = 0.0, T = 10.0)
            r = simular(p, 0.0, 10.0, 401)
            push!(equilibrios, Dict{String, Any}("n" => n, "eta" => η, "max_abs_ell_menos_ell0" => maximum(abs, r.ℓ .- p.ℓ0),
                                                 "max_abs_P" => maximum(abs, r.P)))
        end
    end
    criterios = Dict{String, Any}(
        "monotonia_positividad_cota" => cotas_ok,
        "balance" => all(f["residuo_balance"] < 1e-9 for f in barrido),
        "pendientes_tardias" => all(abs(f["diferencia_ell"]) ≤ 0.05 && abs(f["diferencia_P"]) ≤ 0.05
                                    for f in tendencias if f["T"] == 1e4),
        "equilibrios" => all(f["max_abs_ell_menos_ell0"] < 1e-14 && f["max_abs_P"] < 1e-14 for f in equilibrios))
    Dict{String, Any}("parametros_base" => comodict(base), "T_barrido" => 10.0, "n_pts_barrido" => 401,
                      "eta_tendencias" => 0.0, "U_tendencias" => 1.0, "n_pts_tendencias" => 11, "T_equilibrios" => 10.0,
                      "U_equilibrios" => 0.0, "P0_equilibrios" => 0.0, "n_pts_equilibrios" => 401, "barrido" => barrido,
                      "tendencias" => tendencias, "equilibrios" => equilibrios, "criterios" => criterios)
end

# Largo límite sin inyección. Con U ≡ 0 la fractura se detiene cuando S se anula, en la única raíz ℓ∞ de
# n·η·∫ s^(1-n/2)(c+s²)^(n-1) ds = S0^n (integral de ℓ0 a ℓ∞). Se busca por bisección; la usan la figura y las pruebas.
function ell_limite(par)
    S0 = (par.c + par.ℓ0^2) * par.P0
    g(L) = par.n * par.η * quadgk(s -> s^(1 - par.n / 2) * (par.c + s^2)^(par.n - 1), par.ℓ0, L)[1] - S0^par.n
    lo = hi = par.ℓ0
    k = 0
    # Con η = 0 no hay raíz, y sin tope este bucle no terminaría.
    while g(hi) < 0 && k < 200
        hi *= 1.5
        k += 1
    end
    @assert k < 200 "sin cambio de signo: con este η la fractura no se detiene"
    for _ in 1:200
        m = (lo + hi) / 2
        g(m) < 0 ? (lo = m) : (hi = m)
    end
    ((lo + hi) / 2, S0)
end

# Cota explícita [ℓ0^(3n/2) + (3/2η)S0^n]^(2/(3n)) del largo límite, que es igualdad en n = 1.
cota_ell_limite(par, S0) = (par.ℓ0^(3par.n / 2) + 3 / (2par.η) * S0^par.n)^(2 / (3par.n))

# Datos del plano de fase (ℓ, P) de la figura del objetivo 1, en dos familias. Libre (U ≡ 0), cuatro presiones
# iniciales: cada trayectoria termina en su equilibrio (ℓ∞, 0), y la figura marca ahí la raíz y la cota.
# Forzada (U constante), cuatro caudales: ninguna tiene equilibrio, y el umbral U*, raíz de ℓ(T; U) - ℓd,
# separa los caudales que alcanzan el objetivo de los que no. Existe porque ℓ(T; ·) crece, y eso lo miden las pruebas.
const P0S_OBJ1 = (0.15, 0.3, 0.6, 1.0)
const US_OBJ1 = (0.0, 0.6, 1.2, 2.0)
const TAU_LARGO_OBJ1 = 200.0
const NPTS_OBJ1 = 201

# Malla de la familia libre: uniforme hasta T, donde está el codo del plano de fase, y logarítmica de ahí a
# TAU_LARGO_OBJ1, donde la trayectoria solo se arrastra al equilibrio.
malla_libre_obj1(p) = vcat(range(0.0, p.T, length = 81),
                           exp10.(range(log10(1.2 * p.T), log10(TAU_LARGO_OBJ1), length = 120)))

function datos_figura_obj1()
    p = Par()
    ts = malla_libre_obj1(p)
    libres = [begin
                  q = con(p; P0 = P0)
                  r = simular(q, 0.0, ts)
                  li, S0 = ell_limite(q)
                  (P0 = P0, ℓ = collect(r.ℓ), P = collect(r.P), ell_inf = li, S0 = S0,
                   cota = cota_ell_limite(q, S0), ell_fin = r.ℓ[end], P_fin = r.P[end])
              end for P0 in P0S_OBJ1]
    forzadas = [begin
                    r = simular(p, U, p.T, NPTS_OBJ1)
                    (U = U, ℓ = collect(r.ℓ), P = collect(r.P), ellT = r.ℓ[end], PT = r.P[end])
                end for U in US_OBJ1]
    Ustar = raiz(U -> simular(p, U, p.T, 3).ℓ[end] - p.ℓd, 0.0, p.qmax)
    ru = simular(p, Ustar, p.T, NPTS_OBJ1)
    umbral = (U = Ustar, ℓ = collect(ru.ℓ), P = collect(ru.P), ellT = ru.ℓ[end], PT = ru.P[end])
    (p = p, libres = libres, forzadas = forzadas, umbral = umbral, Ustar = Ustar,
     ellT_Ustar = simular(p, Ustar, p.T, 3).ℓ[end], tau_largo = TAU_LARGO_OBJ1)
end
