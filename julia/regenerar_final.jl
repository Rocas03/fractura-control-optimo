# Regenera los resultados de la ENTREGA FINAL: los json de `final/resultados/`, sus figuras y las cifras
# del informe. `avance/` está entregado (26-oct) y no se toca: lo regenera `avance/regenerar.sh`, y que
# siga dando los mismos bytes es el control de que no rompimos lo entregado.
# Uso, desde la raíz del repositorio: julia/bin/jl --project=julia julia/regenerar_final.jl
using Fractura, JSON, Test

mkpath(RES_FINAL)

o2 = objetivo2();
oac, dac = acuerdo();
odf = diseno_fisher();
# El objetivo 5 (F4) devuelve dos cosas: el resultado para el json y las series del experimento de
# control, que la figura necesita y el json no guarda.
o5, d5 = objetivo5_at2();
for (nombre, o) in (("objetivo 2", o2), ("acuerdo (F2)", oac), ("diseño de Fisher (F3.3)", odf),
                      ("objetivo 5: AT2 (F4)", o5))
    println(nombre, " · criterios:")
    for (k, v) in sort(collect(o["criterios"]); by = first)
        println("  ", rpad(k, 46), v)
    end
end

# La regla 3 del plan: la corrida se detiene si un criterio da falso. Un criterio que falla y se publica
# igual no es un criterio. Los retirados van aparte, con su razón, y no detienen nada.
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

FIGS = joinpath(RAIZ, "final", "figuras")
# La figura del objetivo 1 (F0): el plano de fase con las dos familias de §3.1. No sale de un json
# —el del objetivo 1 es del avance, congelado—, sino de la derivación y del integrador, así que sus
# cifras las escribe el generador y no hay macro que las publique (§9.2 lo declara).
escribir_figura_obj1(datos_figura_obj1(), FIGS);
escribir_figura_acuerdo(dac, FIGS);
escribir_figura_frontera(odf, FIGS);
# La figura del LQ (F1.4) pide las SERIES, que el json no guarda: rehace el experimento con
# `casos_lq`, la misma lista que usa `objetivo2`, para que figura y cifras no se separen.
escribir_figura_lq(datos_figura_lq(), FIGS);
# La figura del objetivo 5 (F4.2 paso 5): el programa de carga óptimo y el campo de daño que
# produce, los dos del MISMO experimento, con las series que `objetivo5_at2` ya calculó.
escribir_figura_at2_control(d5, FIGS);

# Las cifras del informe: un \newcommand por cada número de los cuatro json recién escritos (cifras_informe.jl).
# Es la misma pieza que `cifras.jl` hace para el avance, y la razón de que ninguna cifra del informe se tipee a mano.
escribir_cifras(cifras_informe(), CIFRAS_FINAL; etiqueta = "cifras_informe")
