"""
Fractura hidráulica: control óptimo del caudal de inyección (MA4703, avance del proyecto).

Cada archivo de `src/` es una sección del notebook que se entrega (`notebook/fractura.ipynb`), que `notebook/derivar.py`
arma a partir de estos mismos archivos: no hay otra copia del código. Las rutas de lectura y escritura cuelgan de
`RAIZ`, la raíz del repositorio; el notebook entregado las reemplaza por `salidas/`.

Uso: `julia --project=julia julia/regenerar.jl` (resultados, semilla, cifras y figuras del reporte) y
`julia --project=julia julia/test/runtests.jl` (pruebas).
"""
module Fractura

using OrdinaryDiffEqVerner, NLPModels, NLPModelsIpopt, QuadGK, JSON, StaticArrays, LinearAlgebra, Printf, Random
using PrecompileTools: @setup_workload, @compile_workload

export objetivo1, objetivo3, objetivo4, semilla_lectura, regresion, cifras_reporte, escribir_cifras, datos_figuras,
       escribir_figuras, escribir_datos_beamer, RAIZ, RES,
       nominal_control, nominal_optima, matrices_lineales, fundamental, fundamental_densa, propagador_final,
       gramiano_controlabilidad, gramiano_observabilidad, kalman_congelado, observabilidad_congelada,
       energia_minima, riccati, ganancia_lq, lazo_cerrado, objetivo2, contraste_semilla, nominales_semilla,
       casos_lq, datos_figura_lq, figura_lq, escribir_figura_lq,
       RES_FINAL, cifras_informe, CIFRAS_FINAL, acuerdo, control_pmp, sup_du_directo, NS_ACUERDO, figura_acuerdo, escribir_figura_acuerdo,
       sensibilidades, sens_presion, presion_adim, PARS_SENS, presion_dim, DIM_BASE, PISO_SENS,
       calibre_escala, calibre_vk, mover_log, invariantes, monomio, jacobiano_dim, residuo_nucleo,
       SIGMA_FISHER, fisher, cota_cramer_rao, inflacion, ajustar, experimento_cramer_rao,
       NS_DISENO, TS_DISENO, escalon, fisher_diseno, logdet_fisher, grad_logdet_fisher, psi, grad_psi,
       optimizar_diseno, frontera, Q_PROYECTADO, caja, grad_proyectado, diseno_fisher,
       GAMMAS_FRONTERA, figura_frontera, escribir_figura_frontera,
       ParAT2, paso_malla, nodos_at2, centros_at2, gc_elem, nodos_fijos, alpha_elem, longitud_danada,
       respuesta_u, energia_disipada, hessiano_disipado, qp_alpha, pgs!, qp_cotas!, alternancia_at2,
       resolver_at2, eps_c, alpha_homogenea, sigma_homogenea, AT2_PICO, sigma_pico_homogenea,
       eps_pico_homogenea, grad_S, residuo_at2, hessiano_at2, dresiduo_dd,
       grad_longitud, trabajo_at2, costo_at2, gradiente_at2, optimizar_at2, pesos_objetivo, carga_exacta_at2,
       objetivo5_at2, pico_at2, orden_malla_at2, invariancia_escala_at2, rama_homogenea_at2,
       deficit_terminal_at2,
       gradiente_contra_fd_at2, quiebre_at2, seguimiento_at2, rampa_at2, P_AT2_CONTROL, BLANCOS_AT2,
       NS_MALLA_AT2, TAUS_QUIEBRE_AT2, figura_at2_control, escribir_figura_at2_control,
       ell_limite, cota_ell_limite, datos_figura_obj1, figura_obj1, escribir_figura_obj1,
       P0S_OBJ1, US_OBJ1, TAU_LARGO_OBJ1

const RAIZ = normpath(joinpath(@__DIR__, "..", ".."))

include("modelo.jl")
include("integrador.jl")
include("raices.jl")
include("directo.jl")
include("disparo.jl")
include("comparaciones.jl")
include("hjb.jl")
include("objetivo1.jl")
include("objetivo3.jl")
include("objetivo4.jl")
include("objetivo2.jl")
include("acuerdo.jl")
include("sensibilidad.jl")
include("fisher.jl")
include("diseno.jl")
include("at2.jl")
include("at2_adjunto.jl")
include("objetivo5.jl")
include("semilla.jl")
include("regresion.jl")
include("cifras.jl")
include("cifras_informe.jl")
include("figuras.jl")
include("datos_beamer.jl")

# Precompila la semilla y la E/S (json, regresión, cifras), que no dependen del integrador. Los objetivos 1, 3 y 4 NO
# van aquí: compilados en la precompilación, sus resultados cambian en la última cifra (ulp) respecto del JIT del
# notebook, y con ellos dos cifras del reporte (medido el 14-oct; el mecanismo no está medido). La versión rápida
# tiene que dar exactamente lo mismo que el notebook entregado. Tiempos medidos: sin workload 15,7 s, con este 9,4 s,
# con la corrida completa 3–4 s pero con otras cifras.
@setup_workload begin
    dir = mktempdir()
    @compile_workload begin
        semilla_lectura()
        foreach(n -> (o = JSON.parsefile(joinpath(RES, n * ".json")); fuera_de_tolerancia(o, o)),
                ("objetivo1_dinamica", "objetivo3_pontryagin", "objetivo4_hjb"))
        escribir_cifras(cifras_reporte(RES), joinpath(dir, "cifras.tex"))
    end
end

end
