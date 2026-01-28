# Librairies
using DataFrames
using Distributions
using Roots
using Statistics
using GLM, StatsPlots, StatsBase

function estim_intercept_Y1(df_bool::Vector{Float64}, target::Float64)
    coef = -1
    regression = df_bool * coef

    RR = exp.(regression)
    prob = RR
    prob_reg = clamp.(prob, 0, 1)
    infection = rand.(Binomial.(1, prob_reg))
    prop_estim = mean(infection)
    Intercept_Infection = log(target / prop_estim)
    return (hpv = infection, intercept = Intercept_Infection)
end

function estim_HPV_Y1(df_bool::Vector{Float64}, Intercept_HPV::Float64)

    coef = -1
    regression = df_bool * coef

    # Calcul pour une prévalence 
    prob_reg = exp.(Intercept_HPV .+ regression)

    # Simulation des infections avec le rajout d'un intercept
    HPV = rand.(Binomial.(1,prob_reg))

    return HPV
end


# SimulerT

n = 5000
incidence = [0.055,0.03,0.05,0.025]
souche = ["10","12","20","31"]
incidence_Y1 = 0.22
M = 5000

# Vaccin
modalite_vac = ["Vaccinated","Unvaccinated"]
vec_vac = []

prop_vac = rand(Uniform(50,68))/100
append!(vec_vac, prop_vac)

cols = Symbol.("HPV_" .* souche)
df_Y2 = DataFrame([fill(NaN, n) for _ in cols], cols)
vec_intercept = []

# estimation intercept Y1
for _ in 1:M
    vaccin = wsample(modalite_vac,[prop_vac, 1 - prop_vac],n)
    vaccin = Float64.(vaccin .== "Vaccinated")
    intercept = estim_intercept_Y1(vaccin,incidence_Y1).intercept
    append!(vec_intercept, intercept) 
end

interceptHPV = mean(vec_intercept)

# Y1 est la souche 18
vec_Y1 =[]

biaisDEMA = Float64[]
effet_vaccinDEMA = Float64[]

biaisLOLA = Float64[]
effet_vaccinLOLA = Float64[]

for i in 1:M
    vaccin = wsample(modalite_vac,[prop_vac, 1 - prop_vac],n)
    vaccin = Float64.(vaccin .== "Vaccinated")
     for i in 1:length(souche)
        simulHPV = rand(Bernoulli(incidence[i]), n)
        df_Y2[!, Symbol("HPV_" * souche[i])] = simulHPV
    end 
    Y1 = estim_HPV_Y1(vaccin,interceptHPV)
    append!(vec_Y1, Y1)

     # Etape 4 : Fusion des données
    data_simul  = hcat(DataFrame(Y1=Y1, vaccin=vaccin), df_Y2)

    # print(first(data_simul,5))

    Y2 = sum(Matrix(df_Y2),dims=2)

    data_simul.Y2 = vec(Y2)

    #Méthode par DEMA
    #jsp si on doit laisser ou non l'intercept
    m = lm(@formula(Y1 ~ vaccin + Y2),data_simul) 
    tbl = DataFrame(coeftable(m))

    # colonne nom et coef
    name_col = names(tbl)[1]
    coef_col = names(tbl)[2]
    idx_Y2  = findfirst(==("Y2"), tbl[!, name_col])
    idx_vac = findfirst(==("vaccin"), tbl[!, name_col])
    push!(biaisDEMA,tbl[idx_Y2, coef_col])
    push!(effet_vaccinDEMA,tbl[idx_vac, coef_col])
    
    #Méthode par LOLA
    m1 = glm(@formula(Y1 ~ vaccin),data_simul, Binomial(),LogLink())
    m2 = glm(@formula(Y2 ~ vaccin),data_simul, Poisson())

    tbl1 = DataFrame(coeftable(m1))
    name_col = names(tbl1)[1]
    coef_col = names(tbl1)[2]
    idx_vac = findfirst(==("vaccin"), tbl1[!, name_col])
    effet_vac_biais = tbl1[idx_vac, coef_col]

    tbl2 = DataFrame(coeftable(m2))
    name_col = names(tbl2)[1]
    coef_col = names(tbl2)[2]
    idx_vac = findfirst(==("vaccin"), tbl2[!, name_col])
    biais = tbl2[idx_vac, coef_col]
    push!(biaisLOLA,biais)
    effet_vac = effet_vac_biais - biais
    push!(effet_vaccinLOLA,effet_vac)


end

histogram(effet_vaccinDEMA,
          bins=50,
          xlabel="effetVac DEMA",
          ylabel="Fréquence",
          title="Distribution de l'effet_vac DEMA",
          color=:lightblue)

histogram(effet_vaccinLOLA,
          bins=50,
          xlabel="effet_vac debiaisé ETIEVANT",
          ylabel="Fréquence",
          title="Distribution de l'effet_vac",
          color=:lightblue)