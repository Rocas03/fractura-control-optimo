# F2: el acuerdo entre los tres métodos. ¿Dan lo mismo el directo, el disparo y HJB?
#
# Los tres cómputos ya existían por separado (objetivos 3 y 4) y estaban cruzados de a pares, pero
# nadie los había puesto en los mismos ejes ni medido su convergencia en una sola tabla. Lo que
# aporta este archivo no es un método nuevo: es la verificación cruzada de los tres, que es lo que
# vuelve creíble el óptimo. Los tres caminos son independientes en lo que calculan:
#
#   directo   discretizar y después optimizar: N incógnitas y un IP sobre ellas (Ipopt).
#   disparo   resolver el problema de contorno del PMP: dos incógnitas, z0, y Newton.
#   HJB       programación dinámica en una malla de (ℓ, P, t), y su ley de feedback en lazo cerrado.
#
# Lo que esto prueba y lo que no: un error de MODELO movería a los tres por igual, así que coincidir
# no valida el modelo. Lo que valida es la IMPLEMENTACIÓN, que es donde se cometen los errores que
# nadie ve: un signo en el adjunto, un índice en la malla, una condición de borde. Así se dice.
#
#   acuerdo()          el disparo, la cascada del directo hasta N = 400, los tres niveles de HJB,
#                      la tabla de convergencia, los criterios, y los datos de la figura.
#   control_pmp        U*(t) del extremal como función de t: la referencia de las comparaciones.
#   sup_du_directo     sup|q_k − U*(punto medio del tramo k)|, la discrepancia entre CURVAS.

const NS_ACUERDO = (50, 100, 200, 400)

# Tolerancia del acuerdo en J entre los tres métodos, con la discretización más fina de cada uno.
# NO viene del plan: la fijo acá, el 14-nov, después de la sonda de costo que midió las brechas, y
# lo digo en vez de presentarla como criterio previo. 5e-4 relativo es dos órdenes de magnitud bajo
# la brecha del peor control razonable (la mejor constante está 6,7 % arriba) y un orden sobre la
# brecha del nivel grueso de HJB (2,5e-3), así que distingue «el método converge» de «coincide».
const TOL_ACUERDO_J = 5e-4

# U*(t) del extremal del disparo. `rp["sol"]` es la solución densa del sistema extremal (estado y
# adjunto), y ustar aplica la minimización del Hamiltoniano con la proyección a [0, qmax].
control_pmp(rp, p) = t -> (y = rp["sol"](t); ustar(x2(y), ζ2(y), p))

# El directo entrega un control constante a trozos. Comparado con un control suave el punto justo es
# el PUNTO MEDIO de cada tramo: ahí el error de la constante a trozos es O(h²), el mismo orden que la
# cascada muestra en J. Comparar en los nodos mediría el ancho de la escalera, no el método.
function sup_du_directo(q, p, Ustar)
    t = malla(p, length(q))
    maximum(abs(q[k] - Ustar((t[k] + t[k+1]) / 2)) for k in eachindex(q))
end

baja(v) = all(v[i+1] < v[i] for i in 1:length(v)-1)

function acuerdo(p = Par())
    criterios = Dict{String, Any}()
    rp = resolver(p, SEMILLA, :penalizado)
    Jstar = rp["J"]
    Ustar = control_pmp(rp, p)
    criterios["disparo_converge"] = rp["exito"] && rp["norma_residuo"] < 1e-10

    dir, q = Dict{Int, Any}(), nothing
    for N in NS_ACUERDO
        r = optimizar(p, N; q0 = q === nothing ? zeros(N) : repeat(q, inner = 2))
        q = r.x
        dir[N] = r
    end
    filas_dir = Dict{String, Any}[]
    for (i, N) in enumerate(NS_ACUERDO)
        b = dir[N].fun - Jstar
        ba = i == 1 ? nothing : dir[NS_ACUERDO[i-1]].fun - Jstar
        push!(filas_dir, Dict{String, Any}("N" => N, "J" => dir[N].fun, "brecha" => b,
                                           "razon" => ba === nothing ? nothing : ba / b,
                                           "orden" => ba === nothing ? nothing : log2(ba / b),
                                           "sup_dU" => sup_du_directo(dir[N].x, p, Ustar)))
    end

    filas_hjb, curvas_hjb = Dict{String, Any}[], []
    for (i, (nl, nP, nt, nu)) in enumerate(NIVELES)
        sol = hjb_fractura(p, nl, nP, nt, nu)
        lc = lazo_cerrado(p, sol)
        b = lc.J - Jstar
        ba = i == 1 ? nothing : filas_hjb[i-1]["brecha"]
        push!(filas_hjb, Dict{String, Any}("nl" => nl, "nP" => nP, "nt" => nt, "n_controles" => nu,
                                           "J_lc" => lc.J, "brecha" => b,
                                           "razon" => ba === nothing ? nothing : ba / b,
                                           "orden" => ba === nothing ? nothing : log2(ba / b),
                                           "sup_dU" => maximum(abs.(lc.U .- Ustar.(lc.t)))))
        push!(curvas_hjb, (t = collect(lc.t), U = collect(lc.U)))
    end

    criterios["brechas_directo_bajan"] = baja([f["brecha"] for f in filas_dir])
    criterios["brechas_hjb_bajan"] = baja([f["brecha"] for f in filas_hjb])
    criterios["orden_directo_cerca_de_2"] = all(abs(f["orden"] - 2) ≤ 0.2 for f in filas_dir if f["orden"] !== nothing)
    criterios["sup_dU_directo_baja"] = baja([f["sup_dU"] for f in filas_dir])
    criterios["sup_dU_hjb_baja"] = baja([f["sup_dU"] for f in filas_hjb])
    Js = (Jstar, filas_dir[end]["J"], filas_hjb[end]["J_lc"])
    maxrel = maximum(abs(a - b) / Jstar for a in Js, b in Js)
    criterios["tres_metodos_coinciden_en_J"] = maxrel ≤ TOL_ACUERDO_J

    res = Dict{String, Any}(
        "parametros_base" => comodict(p),
        "J_disparo" => Jstar,
        "disparo" => Dict{String, Any}("exito" => rp["exito"], "norma_residuo" => rp["norma_residuo"],
                                       "H_oscilacion" => rp["H_oscilacion"], "z0" => rp["z0"]),
        "directo" => filas_dir,
        "hjb" => filas_hjb,
        "acuerdo_en_J" => Dict{String, Any}("J_disparo" => Jstar, "J_directo_N$(NS_ACUERDO[end])" => filas_dir[end]["J"],
                                            "J_hjb_lazo_cerrado_fino" => filas_hjb[end]["J_lc"],
                                            "max_diferencia_relativa" => maxrel, "tolerancia" => TOL_ACUERDO_J,
                                            "tolerancia_fijada_en" => "sesión del 14-nov, después de la sonda de costo; no viene del plan"),
        "criterios" => criterios)

    tp = range(0.0, p.T, length = 801)
    tm200 = malla(p, 200)
    med200 = [(tm200[k] + tm200[k+1]) / 2 for k in 1:200]
    datos = (pmp = (t = collect(tp), U = Ustar.(tp)),
             directo = (N = 200, t = collect(tm200), q = collect(dir[200].x)),
             hjb = (nl = NIVELES[end][1], nt = NIVELES[end][3], t = curvas_hjb[end].t, U = curvas_hjb[end].U),
             dif_directo = (t = med200, d = [abs(dir[200].x[k] - Ustar(med200[k])) for k in 1:200]),
             dif_hjb = (t = curvas_hjb[end].t, d = abs.(curvas_hjb[end].U .- Ustar.(curvas_hjb[end].t))))
    res, datos
end
