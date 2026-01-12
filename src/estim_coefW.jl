function estim_coefW(n::Int64)
   W = rand(Uniform(-2, 2), n, 4)
   return W
end