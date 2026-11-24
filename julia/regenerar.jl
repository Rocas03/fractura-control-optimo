# Regenera los resultados del avance: json de los objetivos 1, 3 y 4, semilla, cifras y figuras del reporte.
# Uso, desde la raíz del repositorio: julia --project=julia julia/regenerar.jl
# Cada bloque «# %% <id>» es una celda del notebook de avance; notebook/derivar.py los lee de acá, así que
# esas marcas no se tocan. Las pruebas van aparte, en julia/test/runtests.jl.
using Fractura, JSON, Test

isdir(RES) || error("no está $RES")

# %% c023
o1 = objetivo1();
o1["criterios"]

# %% c024
o3, f3 = objetivo3();
o3["criterios"]

# %% c025
o4, f4 = objetivo4();
o4["criterios"]

# %% c042
sem = semilla_lectura()
open(joinpath(RAIZ, "semilla", "resultados-lectura.json"), "w") do io
    write(io, JSON.json(sem; pretty = 2), "\n")
end
sem["identidades_numericas_bigfloat"]

# %% c036
regresion("objetivo1_dinamica", o1)
regresion("objetivo3_pontryagin", o3)
regresion("objetivo4_hjb", o4);

# %% c037
@testset "criterios" begin
    for (nombre, o) in (("objetivo1", o1), ("objetivo3", o3), ("objetivo4", o4))
        @testset "$nombre" begin
            for (k, v) in sort(collect(o["criterios"]); by = first)
                @testset "$k" begin
                    @test v === true
                end
            end
        end
    end
end;

# %% c044
for (nombre, o) in (("objetivo1_dinamica", o1), ("objetivo3_pontryagin", o3), ("objetivo4_hjb", o4))
    open(joinpath(RES, nombre * ".json"), "w") do io
        write(io, JSON.json(o; pretty = 2), "\n")
    end
end

# %% c046
escribir_cifras(cifras_reporte())

# %% c048
df = datos_figuras(f3, f4, o4)
escribir_figuras(df)
escribir_datos_beamer(df)
