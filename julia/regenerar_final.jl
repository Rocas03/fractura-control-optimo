# Regenera los resultados de la entrega final: los json de `final/resultados/`, sus figuras y las cifras del
# informe. `avance/` es una entrega cerrada y no se toca acá.
# Uso, desde la raíz del repositorio: julia/bin/jl --project=julia julia/regenerar_final.jl
using Fractura, JSON, Test

mkpath(RES_FINAL)

o2 = objetivo2();
oac, dac = acuerdo();
odf = diseno_fisher();
# El objetivo 5 devuelve el resultado para el json y las series del experimento de control, que la figura
# necesita y el json no guarda.
o5, d5 = objetivo5_at2();
for (nombre, o) in (("objetivo 2", o2), ("acuerdo (F2)", oac), ("diseño de Fisher (F3.3)", odf),
                      ("objetivo 5: AT2 (F4)", o5))
    println(nombre, " · criterios:")
    for (k, v) in sort(collect(o["criterios"]); by = first)
        println("  ", rpad(k, 46), v)
    end
end

# La corrida se detiene si un criterio da falso: un criterio que falla y se publica igual no es criterio.
@testset "criterios de la entrega final" begin
    for (nombre, o) in (("objetivo2", o2), ("acuerdo", oac), ("diseno_fisher", odf),
                          ("objetivo5_at2", o5))
        @testset "$nombre" begin
            for (k, v) in sort(collect(o["criterios"]); by = first)
                @testset "$k" begin
                    @test v === true
                end
            end
        end
    end
end;

for (nombre, o) in (("objetivo2_lineal", o2), ("acuerdo_tres_metodos", oac),
                    ("diseno_fisher", odf), ("objetivo5_at2", o5))
    open(joinpath(RES_FINAL, nombre * ".json"), "w") do io
        write(io, JSON.json(o; pretty = 2), "\n")
    end
    println("escrito: ", joinpath("final", "resultados", nombre * ".json"))
end

# El quinto json, la atribución del cruce de Lyapunov (julia/cruce_atribucion.jl). Tiene que reproducir bit
# a bit el error_rel_maximo de objetivo2_lineal.json, así que va después de escribirlo.
include(joinpath(@__DIR__, "cruce_atribucion.jl"))
open(joinpath(RES_FINAL, "cruce_lyapunov_atribucion.json"), "w") do io
    write(io, JSON.json(cruce_atribucion(); pretty = 2), "\n")
end
println("escrito: ", joinpath("final", "resultados", "cruce_lyapunov_atribucion.json"))

FIGS = joinpath(RAIZ, "final", "figuras")
# La figura del objetivo 1 no sale de un json (el del objetivo 1 es del avance, congelado) sino de la
# derivación y del integrador, así que sus cifras las escribe el generador y no hay macro que las publique.
escribir_figura_obj1(datos_figura_obj1(), FIGS);
escribir_figura_acuerdo(dac, FIGS);
escribir_figura_frontera(odf, FIGS);
# La figura del LQ necesita las series, que el json no guarda: rehace el experimento con `casos_lq`, la
# misma lista que usa `objetivo2`.
escribir_figura_lq(datos_figura_lq(), FIGS);
# La figura del objetivo 5: el programa de carga óptimo y el daño que produce, con las series que
# `objetivo5_at2` ya calculó.
escribir_figura_at2_control(d5, FIGS);

# Las cifras del informe: un \newcommand por cada número de los json recién escritos (cifras_informe.jl).
escribir_cifras(cifras_informe(), CIFRAS_FINAL; etiqueta = "cifras_informe")
