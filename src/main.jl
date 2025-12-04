# Fonction pour tirer une proportion a partir d'un IC
using Distributions
function borneSomme(vec_inf::Vector{Float64},vec_sup::Vector{Float64})
    n = length(vec_inf)
    somme = 0
    prop = Vector{Float64}(undef,n)
    while(somme != 1)
        somme = 0
        for i in 1:n
            prop[i] = round(rand(Uniform(vec_inf[i],vec_sup[i])),digits = 2)
            somme += prop[i]
        end
    end
    return prop
end


# Fonction qui généralise la simulation: 

function simul_Var(mod_var::Vector{String}, prop_var::Vector{Float64}, tauxEvolution::Vector{Vector{Float64}}, n::Int64)
    # Population
    modalite_age = ["18-19","20-24","25-34","35-44"]
    prop_age = [0.06,0.18,0.37,0.39]
    simul_age = n * prop_age
    vec = Vector{Float64}(undef,n)
    for i in 1:length(prop_age)
        prob = prop_var * tauxEvolution[i]
        prob = prob / sum(prob)
        vec = push!(vec,wsample(mod_var,prob,simul_age[i]))
    end
    return vec
end
#prop = borneSomme([16.5,70.9,7.2]/100,[20.7,75.5,9.4]/100)

# Fonction sigmoïde
function sigmoid(x::Vector{Float64})
    return ( 1 ./ ( 1 .+ exp.( - x ) ) ) 
end

#sigmoid([12.0,12.0,130])