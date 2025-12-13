module SimulationHPV

using Roots, Statistics, Distributions, DataFrames
export datasets, sigmoid, estim_Confusion, estim_Vaccin, estim_b0, estim_HPV, borneSomme

# Write your package code here.
include("main.jl")
include("estim_b0.jl")
include("estim_Confusion.jl")
include("estim_HPV.jl")
include("estim_proportionIC.jl")
include("estim_Sigmoid.jl")
include("estim_Vaccin.jl")

end
