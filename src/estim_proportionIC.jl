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

#prop = borneSomme([16.5,70.9,7.2]/100,[20.7,75.5,9.4]/100)