# Linealiza el objetivo 2 en dos nominales: el óptimo del directo (N = 50) y el extremal del disparo.
using Fractura, LinearAlgebra, StaticArrays
import Fractura: Par, optimizar, nominal_optima, Nominal, gramiano_controlabilidad, gramiano_observabilidad,
                 cruce_lyapunov_nominal, riccati, lazo_cerrado, casos_lq, resolver, SEMILLA, ustar, x2, ζ2, N_OBJ2

p = Par()
opt = optimizar(p, N_OBJ2)
nd = nominal_optima(p, opt.x)
rp = resolver(p, SEMILLA, :penalizado)
s = rp["sol"]
ns = Nominal(p, τ -> x2(s(τ)), τ -> (y = s(τ); ustar(x2(y), ζ2(y), p)), p.T)

Q = SMatrix{2, 2}(0.0, 0.0, 0.0, p.wp)
Qf = SMatrix{2, 2}(p.wT, 0.0, 0.0, 0.0)

function medir(nom, nombre)
    λW = eigvals(Symmetric(Matrix(gramiano_controlabilidad(nom))))
    λM = eigvals(Symmetric(Matrix(gramiano_observabilidad(nom))))
    cr = cruce_lyapunov_nominal(nom)
    P = riccati(nom; Q, R = p.wq, Qf)
    ℓT̄ = nom.x(p.T)[1]
    fac = Dict{String, Float64}()
    for (cual, x0, pp) in casos_lq(p)
        ab = lazo_cerrado(nom, P; x0, p = pp, R = p.wq, lazo_abierto = true)
        ce = lazo_cerrado(nom, P; x0, p = pp, R = p.wq)
        fac[cual] = abs(ab.xT[1] - ℓT̄) / abs(ce.xT[1] - ℓT̄)
    end
    println(nombre)
    println("  ell_T = ", ℓT̄, "   P_T = ", nom.x(p.T)[2])
    println("  lambda_W = ", λW, "   cond_W = ", λW[2] / λW[1])
    println("  lambda_M = ", λM, "   cond_M = ", λM[2] / λM[1])
    println("  cruce Lyapunov error_rel_max = ", cr.error_rel_maximo)
    println("  factor LQ: P0 = ", fac["P0"], "   eta = ", fac["eta"])
    (λW, λM, fac)
end

println("J directo = ", opt.fun, "   J disparo = ", rp["J"], "   brecha = ", opt.fun - rp["J"],
        "   disparo exito = ", rp["exito"])
a = medir(nd, "NOMINAL = directo N=50")
b = medir(ns, "NOMINAL = extremal del disparo")
println("diferencia relativa lambda_W: ", abs.(a[1] .- b[1]) ./ abs.(b[1]))
println("diferencia relativa lambda_M: ", abs.(a[2] .- b[2]) ./ abs.(b[2]))
println("diferencia relativa factor: P0 ", abs(a[3]["P0"] - b[3]["P0"]) / b[3]["P0"],
        "  eta ", abs(a[3]["eta"] - b[3]["eta"]) / b[3]["eta"])
