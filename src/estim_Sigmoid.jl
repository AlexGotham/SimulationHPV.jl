# Fonction sigmoïde
using Distributions
function sigmoid(x::Float64)
    return ( 1 ./ ( 1 .+ exp.( - x ) ) ) 
end

