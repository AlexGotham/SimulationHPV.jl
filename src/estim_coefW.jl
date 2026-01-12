function estim_coefW(n::Int64)
    W = matrix{Float64}(undef, n, 4)
    for i in 1:n
        W[i,1] = rand(Uniform(-2,2))
        W[i,2] = rand(Uniform(-2,2))
        W[i,3] = rand(Uniform(-2,2))
        W[i,4] = rand(Uniform(-2,2))
    end
    return W
end