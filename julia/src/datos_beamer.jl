# Curvas del Beamer animado del avance (avance/video/beamer/datos/*.csv), desde los mismos vectores que dibujan las
# figuras del reporte y con el mismo formato %.6g, así que cada CSV repite las coordenadas de su figura. El Beamer los
# lee con pgfplotstable (comun.tex). Reemplaza a avance/video/beamer/extraer.py, que los sacaba de las figuras.
const DATOS_BEAMER = joinpath(RAIZ, "avance", "video", "beamer", "datos")

function csv_beamer(archivo, cabecera, columnas...)
    open(archivo, "w") do io
        println(io, cabecera)
        for fila in zip(columnas...)
            println(io, join((@sprintf("%.6g", x) for x in fila), ","))
        end
    end
end

function escribir_datos_beamer(d, dir = DATOS_BEAMER)
    mkpath(dir)
    nombres = ("extremal", "mejorconst", "igualvol", "rampa")
    length(d.comparacion) == length(nombres) || error("datos_beamer: se esperaban 4 curvas de comparación")
    for (nom, c) in zip(nombres, d.comparacion)
        csv_beamer(joinpath(dir, nom * ".csv"), "t,U,ell,P", c[2], c[3], c[4], c[5])
    end
    csv_beamer(joinpath(dir, "cota-disparo.csv"), "t,U", d.cota.t, d.cota.U)
    csv_beamer(joinpath(dir, "cota-directo.csv"), "t,U", d.cota.td, d.cota.Ud)
    csv_beamer(joinpath(dir, "retro-pmp.csv"), "t,U", d.retro.t, d.retro.U)
    for n in d.retro.niveles
        csv_beamer(joinpath(dir, "retro-hjb$(n.nl).csv"), "t,U", n.t, n.U)
    end
    println("datos del Beamer: $(length(readdir(dir))) CSV en $(basename(dir))/")
end
