# Fonction Infection pour toutes souches et vaccin

using DataFrames

function estim_HPV(df_bool::DataFrame, W::Vector{Float64}, souche::Int64, target::Float64 )
    A = [- 0.5, 1, 1.5, 2, -0.5, 1.5, 2.5, -0.5, 1, 1.5, 2, 2.5] # coeff de confusion fixé
    cible = [6, 11, 16,18]

    # Vaccin prend -10 si souche concernée, sinon 0
    if souche in cible
        coefVac = -10
    else
        coefVac = 0
    end
    coef = vcat(coefVac, W, A)
    regression = Matrix(df_bool) * coef

    # Calcul pour une prévalence 
    # β₀ = estim_b0.(target,regression)
    β₀ = estim_b0_log.(target,regression)
    # prob_reg = sigmoid.(regression+β₀)
    prob_reg = exp.(regression+β₀)

    # Simulation des infections avec le rajout d'un intercept
    HPV = rand.(Binomial.(1,prob_reg))

    return HPV
end


