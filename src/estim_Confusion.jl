# Fonction qui généralise la simulation: 
using Distributions
function Estim_Confusion(mod_var::Vector{String}, prop_var::Vector{Float64}, tauxEvolution::Vector{Vector{Float64}}, n::Int64, proportion_age::Vector{Float64})
    # Population
    modalite_age = ["18-19","20-24","25-34","35-44"]
    prop_age = proportion_age
    simul_age = round.(Int,n * prop_age)
    vec = Vector{String}(undef,0)
    for i in 1:length(prop_age)
        prob = prop_var .* tauxEvolution[i]
        prob = prob / sum(prob)
        vec = vcat(vec,wsample(mod_var,prob,simul_age[i]))
    end
    return vec
end

# taux_age_sex = [[0.01,0.03 , 0.02, 0],[1,1.6,1,0.4],[1,1.4,1,0.6],[1,1.4,1,0.6]]
# propo_age = [0.06,0.18,0.37,0.39]
# moda_ageSexe = ["18-19","20-24","25-34","35-44"]
# propo_ageSexe = [0.1,0.2,0.4,0.3]
# n = 5
# Estim_Confusion(moda_ageSexe,propo_ageSexe,taux_age_sex,n,propo_age)