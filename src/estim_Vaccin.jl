# Fonction Infection pour toutes souches et vaccin

using DataFrames
using Distributions

function estim_Vaccin(df_bool::DataFrame)
    A = [- 0.5, 1, 1.5, 2, -0.5, 1.5, 2.5, -0.5, 1, 1.5, 2, 2.5] # coeff de confusion fixé
    regression = Matrix(df_bool) * A

    # Calcul pour une prévalence 
    # β₀ = estim_b0.(target,regression)
    prob_reg = sigmoid.(regression)

    # Simulation des infections avec le rajout d'un intercept
    Vaccin = rand.(Binomial.(1,prob_reg))

    return Vaccin
end


df_bool = DataFrame(rand(Bool, 100, 12),:auto)
print(estim_Vaccin(df_bool))