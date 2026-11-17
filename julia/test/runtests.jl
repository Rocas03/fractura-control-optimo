# Pruebas del paquete Fractura: julia --project=julia julia/test/runtests.jl (o Pkg.test() con julia/ activo).
# ayudas.jl y pruebas.jl son, tal cual, las dos celdas «Pruebas» del notebook que se entrega, donde todo el código está
# en el mismo espacio de nombres; aquí se importan los nombres del paquete, exportados o no, para que corran igual.
using Fractura
using OrdinaryDiffEqVerner, NLPModels, NLPModelsIpopt, QuadGK, JSON, StaticArrays, LinearAlgebra, Printf, Random, Test
for n in names(Fractura; all = true)
    s = string(n)
    (startswith(s, "#") || n in (:Fractura, :eval, :include)) && continue
    @eval import Fractura: $n
end

include("ayudas.jl")
include("pruebas.jl")
