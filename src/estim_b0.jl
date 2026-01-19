using Roots
using Statistics

# Fonction pour estimer beta0 en fonction de l'incidence

function estim_b0(incidence::Float64, reg::Float64)
   function β₀(β) 
    mean(sigmoid.(reg.+β)) - incidence
   end
   res = find_zero(β₀,[-100,100])
   return res
end

function estim_b0_log(incidence::Float64, reg::Float64)
   function β₀(β) 
    mean(exp.(reg .+ β)) - incidence
   end
   res = find_zero(β₀, [-20.0,0.0])
   return res
end


