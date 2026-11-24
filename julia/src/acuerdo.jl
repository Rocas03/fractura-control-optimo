# Acuerdo de los tres métodos en una misma tabla de convergencia: el directo (discretizar y optimizar con Ipopt),
# el disparo (Newton sobre ζ(0) en el sistema de Pontryagin) y el lazo cerrado de HJB (programación dinámica).
# Calculan por caminos independientes, así que coincidir valida la implementación (un signo en el adjunto, un
# índice en la malla, un borde), pero no el modelo: un error de modelo los movería a los tres por igual.

const NS_ACUERDO = (50, 100, 200, 400)

# Tolerancia relativa del acuerdo en J, con la discretización más fina de cada método. Se fijó después de medir
# las brechas, no antes. Queda dos órdenes bajo la brecha de la mejor constante (6,7 %) y bajo la del nivel
# grueso de HJB (2,5e-3), así que separa coincidir de solo parecerse.
const TOL_ACUERDO_J = 5e-4

# U*(t) del extremal del disparo, la referencia de las comparaciones: se evalúan estado y adjunto en t y se
# aplica U* = proy(ζP / (wq D)).
control_pmp(rp, p) = t -> (y = rp["sol"](t); ustar(x2(y), ζ2(y), p))

# sup |q_k - U*| con U* evaluado en el punto medio de cada tramo: ahí el error de la constante a trozos es O(h²),
# el mismo orden que muestra J; en los nodos se mediría el ancho de la escalera y no el método.
function sup_du_directo(q, p, Ustar)
    t = malla(p, length(q))
    maximum(abs(q[k] - Ustar((t[k] + t[k+1]) / 2)) for k in eachindex(q))
end

baja(v) = all(v[i+1] < v[i] for i in 1:length(v)-1)

# Disparo, cascada del directo hasta N = 400 y los tres niveles de HJB, con brechas en J, orden observado,
# criterios y los datos de la figura.
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
