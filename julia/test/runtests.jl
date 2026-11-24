# Pruebas del paquete: julia --project=julia julia/test/runtests.jl, o Pkg.test() con julia/ activo.
# ayudas.jl y pruebas.jl son también las celdas de pruebas del notebook, donde todo vive en un mismo espacio de
# nombres; por eso aquí se importan todos los nombres del paquete, exportados o no.
using Fractura
using OrdinaryDiffEqVerner, NLPModels, NLPModelsIpopt, QuadGK, JSON, StaticArrays, LinearAlgebra, Printf, Random, Test
for n in names(Fractura; all = true)
    s = string(n)
    (startswith(s, "#") || n in (:Fractura, :eval, :include)) && continue
    @eval import Fractura: $n
end

include("ayudas.jl")
include("pruebas.jl")
