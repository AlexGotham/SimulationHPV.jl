# Fonction Infection pour toutes souches et vaccin

using DataFrames
using Distributions

function estim_Vaccin(df_bool::DataFrame, prop_vaccin::Float64)
    A = [- 0.5, 1, 1.5, 2, -0.5, 1.5, 2.5, -0.5, 1, 1.5, 2, 2.5] # coeff de confusion fixé
    regression = Matrix(df_bool) * A

    # Calcul pour une prévalence 
    β₀ = estim_b0.(prop_vaccin,regression)
    prob_reg = sigmoid.(regression + β₀)

    # Simulation des infections avec le rajout d'un intercept
    Vaccin = rand.(Binomial.(1,prob_reg))

    return Vaccin
end


# df_bool_vaccin = DataFrame(rand(Bool, 100, 12),:auto)
# prop_vaccin_simul = 0.55
# print(estim_Vaccin(df_bool, prop_vaccin_simul))