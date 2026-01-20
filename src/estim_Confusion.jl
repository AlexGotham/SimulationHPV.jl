# Fonction qui généralise la simulation: 
using Distributions
function estim_Confusion(mod_var::Vector{String}, prop_var::Vector{Float64}, tauxEvolution::Vector{Vector{Float64}}, n::Int64, proportion_age::Vector{Float64})
    # Population
    # modalite_age = ["18-19","20-24","25-34","35-44"]
    # prop_age = proportion_age
    # simul_age = round.(Int,n * prop_age)
    simul_age = rand(Multinomial(n, proportion_age))
    vec = Vector{String}(undef,0)
    for i in 1:length(proportion_age)
        prob = prop_var .* tauxEvolution[i]
        prob = prob / sum(prob)
        # vec = vcat(vec,wsample(mod_var,prob,simul_age[i]))
        append!(vec, wsample(mod_var, prob, simul_age[i]))
    end
    return vec
end

n = 1000
modalite_age = ["18-19","20-24","25-34","35-44"]
born_inf_age = [5,16,34.7,37.1]/100
born_sup_age = [6.8,19.2,39.2, 42.3]/100
# prop_age = borneSomme(born_inf_age,born_sup_age)

# Avec Dirichlet
k = 80
moy_age = (born_inf_age .+ born_sup_age) ./2
alpha_age = moy_age .* k
prop_age = estim_dirichlet(alpha_age,born_inf_age, born_sup_age)

proportion_age = prop_age
simul_age = rand(Multinomial(n, proportion_age))
sum(simul_age)  