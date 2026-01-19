using Distributions
function estim_dirichlet(alpha, born_inf, born_sup)
    dist = Dirichlet(alpha)
    while true
        p = rand(dist)
        if all(born_inf .<= p .<= born_sup)
            return p
        end
    end
end