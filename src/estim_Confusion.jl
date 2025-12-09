# Fonction qui généralise la simulation: 
using Distributions
function estim_Confusion(mod_var::Vector{String}, prop_var::Vector{Float64}, tauxEvolution::Vector{Vector{Float64}}, n::Int64, proportion_age::Vector{Float64})
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

