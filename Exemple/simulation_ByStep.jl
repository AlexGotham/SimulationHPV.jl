# Librairies
using DataFrames
using Distributions
using Roots
using Statistics
using GLM, StatsPlots, StatsBase

# Présentations du code de ce matin : 
# - Fonctions utilisées pour tout le code
# - Step 1 : Simulation des proportions des variables
# - Step 2 : Simulation de l'intercept du Vaccin (car effet de confusion de W et A sur T)
# - Step 3 : Simulation de l'intercept des Infections
# - Step 4 : Utilisation des simulations précédentes pour simuler les méthodes

# ---------------------------------------------------
################### FONCTIONS #######################
# ---------------------------------------------------
# Fonction d'estimation de l'intercept des Infections pour toutes souches et si l'individu est vacciné ou non vaccin
function estim_intercept_HPV(df_bool::DataFrame, W::Vector{Float64}, souche::Int64, target::Float64)
    A = [- 0.5, 1, 1.5, 2, -0.5, 1.5, 2.5, -0.5, 1, 1.5, 2, 2.5] # coeff de confusion fixé
    cible = [6, 11, 16,18]
    # Vaccin prend Intercept_Vac si souche concernée, sinon 0
    if souche in cible
        coefVac = -100
    else
        coefVac = 0
    end
    coef = vcat(coefVac, W, A) ./100
    regression = Matrix(df_bool) * coef
    # Estimation de l'intercept de l'infection
    RR = exp.(regression)
    prob = RR
    prob_reg = clamp.(prob, 0, 1)
    infection = rand.(Binomial.(1, prob_reg))
    prop_estim = mean(infection)
    Intercept_Infection = log(target / prop_estim)
    return (hpv = infection, intercept = Intercept_Infection)
end

# Fonction d'application des intercepts des Infections déterminés précédemment
function estim_HPV_log(df_bool::DataFrame, W::Vector{Float64}, souche::Int64, Intercept_HPV::Float64)
    A = [- 0.5, 1, 1.5, 2, -0.5, 1.5, 2.5, -0.5, 1, 1.5, 2, 2.5] # coeff de confusion fixé
    cible = [6, 11, 16,18]
    # Vaccin prend Intercept_Vac si souche concernée, sinon 0
    if souche in cible
        coefVac = -100
    else
        coefVac = 0
    end
    coef = vcat(coefVac, W, A) ./100
    regression = Matrix(df_bool) * coef
    # Calcul pour une prévalence selon l'intercept
    prob_reg = exp.(Intercept_HPV .+ regression)
    # Simulation des infections avec le rajout d'un intercept
    HPV = rand.(Binomial.(1,prob_reg))
    return HPV
end

# Fonction pour tirer une proportion a partir d'un IC Par Dirichlet
function estim_dirichlet(alpha, born_inf, born_sup)
    dist = Dirichlet(alpha)
    while true
        p = rand(dist)
        # Vérification de la présence dans l'intervalle de confiance
        if all((born_inf .<= p) .& (p .<= born_sup))
            return p
        end
    end
end

# Fonction de génération de Coef W aléatoire
function estim_coefW(n::Int64)
   W = rand(Uniform(-2, 2), n, 4)
   return W
end

# Fonction qui généralise la simulation des variables avec un effet de confusion de l'âge: 
function estim_Confusion(mod_var::Vector{String}, prop_var::Vector{Float64}, tauxEvolution::Vector{Vector{Float64}}, n::Int64, proportion_age::Vector{Float64})
    simul_age = rand(Multinomial(n, proportion_age))
    vec = Vector{String}(undef,0)
    for i in 1:length(proportion_age)
        prob = prop_var .* tauxEvolution[i]
        prob = prob / sum(prob)
        append!(vec, wsample(mod_var, prob, simul_age[i]))
    end
    return vec
end

# Fonction d'estimation de l'intercept du vaccin avec les confusions de A et W
function estim_Intercept_Vaccin_log(df_bool::DataFrame, prop_vaccin::Float64)
    W = [2.5,-0.5,1,2.5] # effet de l'age sur le vaccin
    A = [- 0.5, 1, 1.5, 2, -0.5, 1.5, 2.5, -0.5, 1, 1.5, 2, 2.5] # coeff de confusion fixé
    coef = vcat(W,A) ./100
    # Calcul pour une prévalence lié à la prop
    Regression = Matrix(df_bool) * coef
    RR = exp.(Regression)
    prob = RR
    prob_reg = clamp.(prob, 0, 1)
    Vaccin = rand.(Binomial.(1, prob_reg))
    prop_estim = mean(Vaccin)
    Intercept_T = log(prop_vaccin / prop_estim)
    return (intercept = Intercept_T , vaccin = Vaccin)
end

# Fonction d'application de l'intercept du Vaccin pour la simulation de la variable de Vaccination
function estim_Vaccin_log(df_bool::DataFrame, Intercept_T::Float64)
    W = [2.5,-0.5,1,2.5] # effet de l'age sur le vaccin
    A = [- 0.5, 1, 1.5, 2, -0.5, 1.5, 2.5, -0.5, 1, 1.5, 2, 2.5] # coeff de confusion fixé
    coef = vcat(W,A) ./ 100
    regression = Matrix(df_bool)* coef
    # Calcul pour une prévalence 
    prob_reg = exp.(Intercept_T .+ regression)
    # Simulation des infections avec le rajout d'un intercept
    Vaccin = rand.(Binomial.(1,prob_reg))
    return Vaccin
end

# ---------------------------------------------------
##################### STEP 1 ########################
# ---------------------------------------------------
# Estimation des proportions sur 5000 itérations
n = 1000
target = [0.22,0.055,0.03,0.05,0.025]
souche = ["18","10","12","20","31"]
M = 5000

# Initialisation des stocks
# Age
modalite_age = ["18-19","20-24","25-34","35-44"]
born_inf_age = [5,16,34.7,37.1]/100
born_sup_age = [6.8,19.2,39.2, 42.3]/100
df_age = DataFrame(
    Symbol.(modalite_age) .=> eachcol(zeros(M, 4))
)
# Vaccin
modalite_vac = ["Vaccinated","Unvaccinated"]
vec_vac = []
# Age first sex
modalite_age_sex = ["13-15","16-17","18-19","20 et +"]
born_inf_age_sex = [24.1,40.7,17.3,9.6]/100
born_sup_age_sex = [28.1,45.8,21.6,13.3]/100
df_age_sex = DataFrame(
    Symbol.(modalite_age_sex) .=> eachcol(zeros(M, 4))
)
# Condomless
modalite_condomless = ["0","1","2 et +"]
born_inf_condomless = [16.5,70.9,7.2]/100
born_sup_condomless = [20.7,75.5,9.4]/100
df_condomless = DataFrame(
    Symbol.(modalite_condomless) .=> eachcol(zeros(M, 3))
)
# Partner
modalite_partner = ["0","1","2","3-4","5 et +"]
born_inf_partner = [0.9,54.9,12.7,11.4,12.8]/100
born_sup_partner = [2.4,59.8,15.8,14.4,15.6]/100
df_partner = DataFrame(
    Symbol.(modalite_partner) .=> eachcol(zeros(M, 5))
)

# Boucle Population
for i in 1:M
    # Pour Dirichlet
    k = 80
    # ----------------------------------------------------------------------------
    # Etape 1 : Simulation des variables
    # simulation age Avec Dirichlet
    moy_age = (born_inf_age .+ born_sup_age) ./2
    alpha_age = moy_age .* k
    df_age[i, :] = estim_dirichlet(alpha_age, born_inf_age, born_sup_age)

    # # simulation  HPV vaccination status
    prop_vac_i = rand(Uniform(50,68))/100
    append!(vec_vac, prop_vac_i)

    # simulation age first sex
    moy_age_sex = (born_inf_age_sex .+ born_sup_age_sex) ./2
    alpha_age_sex = moy_age_sex .* k
    df_age_sex[i, :] = estim_dirichlet(alpha_age_sex,born_inf_age_sex, born_sup_age_sex)

    # proportion condomless sex with partner, past year
    moy_condomless = (born_inf_condomless .+ born_sup_condomless) ./2
    alpha_condomless = moy_condomless .* k
    df_condomless[i, :] = estim_dirichlet(alpha_condomless,born_inf_condomless, born_sup_condomless)

    # proportion de Total partners, past 5 years
    moy_partner = (born_inf_partner .+ born_sup_partner) ./2
    alpha_partner = moy_partner.* k
    df_partner[i, :] = estim_dirichlet(alpha_partner,born_inf_partner, born_sup_partner)  
end

PROP_AGE = mean.(eachcol(df_age))
PROP_VAC = mean(vec_vac)
PROP_AGE_SEX = mean.(eachcol(df_age_sex))
PROP_CONDOMLESS = mean.(eachcol(df_condomless))
PROP_PARTNER = mean.(eachcol(df_partner))

println("Verif Age : ", sum(PROP_AGE))
println("Verif Age Sex : ", sum(PROP_AGE_SEX))
println("Verif Condomless : ", sum(PROP_CONDOMLESS))
println("Verif Partner : ", sum(PROP_PARTNER))


# ---------------------------------------------------
################### RESULTATS #######################
# ---------------------------------------------------
# Sauvegarde des résultats en brut !!!!!!!!!!!!!!!!!!
# Age
proportion_age = [0.05845492872805624,
0.17571140375916522,
0.3694368528694803,
0.39639681464329823]
"Verif Age : 1.0"

# Age Sexe
proportion_age_sex = [0.2606516856419951,
0.43219808926791065,
0.1933099374922233,
0.11384028759787101]
"Verif Age Sex : 1.0"

# Condomless
proportion_condomless = [0.18516910637152836,
0.7323120679567233,
0.08251882567174829]
"Verif Condomless : 1.0"

# Partner
proportion_partner = [0.015330633164757393,
0.5729406526715017,
0.14189934961754136,
0.12838568114493593,
0.14144368340126356]
"Verif Partner : 1.0"

PROP_VAC = 0.5906981405248859


# ---------------------------------------------------
##################### STEP 2 ########################
# ---------------------------------------------------
# Estimation des intercept sur 5000 itérations avec les proportions validées
# ----------------------------------------------------------------------------
intercept_vaccin = []
for i in 1:M
    # Etape 2 : Simulation des variables avec interaction de l'âge
    # Variable Âge
    age = wsample(modalite_age, proportion_age, n) 

    # Première fois en fonction de l'âge
    taux_age_sex = [
        [1.2, 1.6, 1.2, 0], 
        [1,1.6,1,0.4],  
        [1,1.4,1,0.6],   
        [1,1.4,1,0.6]
    ]
    SimulFirstTime = estim_Confusion(modalite_age_sex, proportion_age_sex, taux_age_sex , n, proportion_age)

    # Utilisation de protection en fonction de l'âge
    taux_condomless = [
        [1.2, 1.6, 1.2], 
        [1,1.6,1],  
        [1,1.4,1],   
        [1,1.4,1]
    ]
    SimulCondom = estim_Confusion(modalite_condomless, proportion_condomless, taux_condomless, n, proportion_age)

    # Nombre de partenaire en fonction de l'âge
    taux_partner = [
        [0, 2, 1.2,1,0.8], 
        [0.3,1.8,1,1,0.9],  
        [0.3,1.5,1.1,1.1,1],   
        [0.3,1.7,1.1,1.1,0.8]
    ]
    SimulPartner = estim_Confusion(modalite_partner, proportion_partner, taux_partner, n, proportion_age)

    # Regroupement des données dans un dataframe
    # Age = sort(age)
    DF = DataFrame(Age = age,  SimulFirstTime = SimulFirstTime, SimulCondom = SimulCondom, SimulPartner = SimulPartner)

    # ----------------------------------------------------------------------------
    # Etape 3 : Infection pour n souches
    # Binariser toutes les variables
    DF_bool = DataFrame(
    age_18_19 = ifelse.(DF.Age .== "18-19",1,0),
    age_20_24 = ifelse.(DF.Age .== "20-24",1,0),
    age_25_34 = ifelse.(DF.Age .== "25-34",1,0),
    age_35_44 = ifelse.(DF.Age .== "35-44",1,0),
    first_time_13_15 = ifelse.(DF.SimulFirstTime .== "13-15",1,0),
    first_time_16_17 = ifelse.(DF.SimulFirstTime .== "16-17",1,0),
    first_time_18_19 = ifelse.(DF.SimulFirstTime .== "18-19",1,0),
    first_time_20 = ifelse.(DF.SimulFirstTime .== "20 et +",1,0),
    condomless_0 = ifelse.(DF.SimulCondom .== "0",1,0),
    condomless_1 = ifelse.(DF.SimulCondom .== "1",1,0),
    condomless_2 = ifelse.(DF.SimulCondom .== "2 et +",1,0),
    partner_0 = ifelse.(DF.SimulPartner .== "0",1,0),
    partner_1 = ifelse.(DF.SimulPartner .== "1",1,0),
    partner_2 = ifelse.(DF.SimulPartner .== "2",1,0),
    partner_3_4 = ifelse.(DF.SimulPartner .== "3-4",1,0),
    partner_5 = ifelse.(DF.SimulPartner .== "5 et +",1,0)  
    )

    # Simulation des vaccins en fonction des covariables
    res = estim_Intercept_Vaccin_log(DF_bool, PROP_VAC)
    intercept_vaccin_i = res.intercept
    append!(intercept_vaccin, intercept_vaccin_i)
end

# ---------------------------------------------------
################### RESULTATS #######################
# ---------------------------------------------------
# Sauvegarde des résultats en brut !!!!!!!!!!!!!!!!!!
Intercept_T = mean(intercept_vaccin)



# ---------------------------------------------------
##################### STEP 3 ########################
# ---------------------------------------------------
# Simulation des infections (FONCTION)

cols = Symbol.("HPV_" .* souche)
df_infection = DataFrame([fill(NaN, M) for _ in cols], cols)

for i in 1:M
    # Variable Âge
    age = wsample(modalite_age, proportion_age, n) 

    # Première fois en fonction de l'âge
    taux_age_sex = [
        [1.2, 1.6, 1.2, 0], 
        [1,1.6,1,0.4],  
        [1,1.4,1,0.6],   
        [1,1.4,1,0.6]
    ]
    SimulFirstTime = estim_Confusion(modalite_age_sex, proportion_age_sex, taux_age_sex , n, proportion_age)

    # Utilisation de protection en fonction de l'âge
    taux_condomless = [
        [1.2, 1.6, 1.2], 
        [1,1.6,1],  
        [1,1.4,1],   
        [1,1.4,1]
    ]
    SimulCondom = estim_Confusion(modalite_condomless, proportion_condomless, taux_condomless, n, proportion_age)

    # Nombre de partenaire en fonction de l'âge
    taux_partner = [
        [0, 2, 1.2,1,0.8], 
        [0.3,1.8,1,1,0.9],  
        [0.3,1.5,1.1,1.1,1],   
        [0.3,1.7,1.1,1.1,0.8]
    ]
    SimulPartner = estim_Confusion(modalite_partner, proportion_partner, taux_partner, n, proportion_age)

    # Regroupement des données dans un dataframe
    # Age = sort(age)
    DF = DataFrame(Age = age,  SimulFirstTime = SimulFirstTime, SimulCondom = SimulCondom, SimulPartner = SimulPartner)

    # ----------------------------------------------------------------------------
    # Etape 3 : Infection pour n souches
    # Binariser toutes les variables
    DF_bool = DataFrame(
    age_18_19 = ifelse.(DF.Age .== "18-19",1,0),
    age_20_24 = ifelse.(DF.Age .== "20-24",1,0),
    age_25_34 = ifelse.(DF.Age .== "25-34",1,0),
    age_35_44 = ifelse.(DF.Age .== "35-44",1,0),
    first_time_13_15 = ifelse.(DF.SimulFirstTime .== "13-15",1,0),
    first_time_16_17 = ifelse.(DF.SimulFirstTime .== "16-17",1,0),
    first_time_18_19 = ifelse.(DF.SimulFirstTime .== "18-19",1,0),
    first_time_20 = ifelse.(DF.SimulFirstTime .== "20 et +",1,0),
    condomless_0 = ifelse.(DF.SimulCondom .== "0",1,0),
    condomless_1 = ifelse.(DF.SimulCondom .== "1",1,0),
    condomless_2 = ifelse.(DF.SimulCondom .== "2 et +",1,0),
    partner_0 = ifelse.(DF.SimulPartner .== "0",1,0),
    partner_1 = ifelse.(DF.SimulPartner .== "1",1,0),
    partner_2 = ifelse.(DF.SimulPartner .== "2",1,0),
    partner_3_4 = ifelse.(DF.SimulPartner .== "3-4",1,0),
    partner_5 = ifelse.(DF.SimulPartner .== "5 et +",1,0)  
    )

    # Simulation des vaccins en fonction des covariables
    SimulVaccin = estim_Vaccin_log(DF_bool, Intercept_T)
    DF.SimulVaccin = SimulVaccin
    DF_bool[!, :Vaccin] = SimulVaccin

    # Paramètres de l'âge
    W = estim_coefW(length(souche))

    # Génération des infections par souches
    for j in 1 : length(souche)
        infection = estim_intercept_HPV(DF_bool, W[j,:], parse(Int64, souche[j]), target[j])
        Intercept_infection = infection.intercept
        df_infection[i, "HPV_" * souche[j]] = Intercept_infection
    end 
end
# ---------------------------------------------------
################### RESULTATS #######################
# ---------------------------------------------------
# Sauvegarde des résultats en brut !!!!!!!!!!!!!!!!!!
intercept_infection = mean.(eachcol(df_infection))


# ---------------------------------------------------
##################### STEP 4 ########################
# ---------------------------------------------------
M = 10000
biaisDEMA = Float64[]
effet_vaccinDEMA = Float64[]
lowerDEMA = Float64[]
upperDEMA = Float64[]
nbDEMA = 0

biaisLOLA = Float64[]
effet_vaccinLOLA = Float64[]
lowerLOLA = Float64[]
upperLOLA = Float64[]
nbLOLA = 0

for i in 1:M
    # Simulation des jeux de données pour chaque méthodes
    cols = Symbol.("HPV_" .* souche)
    df_HPV = DataFrame([fill(NaN, n) for _ in cols], cols)
    # Variable Âge
    age = wsample(modalite_age, proportion_age, n) 

    # Première fois en fonction de l'âge
    taux_age_sex = [
        [1.2, 1.6, 1.2, 0], 
        [1,1.6,1,0.4],  
        [1,1.4,1,0.6],   
        [1,1.4,1,0.6]
    ]
    SimulFirstTime = estim_Confusion(modalite_age_sex, proportion_age_sex, taux_age_sex , n, proportion_age)

    # Utilisation de protection en fonction de l'âge
    taux_condomless = [
        [1.2, 1.6, 1.2], 
        [1,1.6,1],  
        [1,1.4,1],   
        [1,1.4,1]
    ]
    SimulCondom = estim_Confusion(modalite_condomless, proportion_condomless, taux_condomless, n, proportion_age)

    # Nombre de partenaire en fonction de l'âge
    taux_partner = [
        [0, 2, 1.2,1,0.8], 
        [0.3,1.8,1,1,0.9],  
        [0.3,1.5,1.1,1.1,1],   
        [0.3,1.7,1.1,1.1,0.8]
    ]
    SimulPartner = estim_Confusion(modalite_partner, proportion_partner, taux_partner, n, proportion_age)

    # Regroupement des données dans un dataframe
    # Age = sort(age)
    DF = DataFrame(Age = age,  SimulFirstTime = SimulFirstTime, SimulCondom = SimulCondom, SimulPartner = SimulPartner)

    # ----------------------------------------------------------------------------
    # Etape 3 : Infection pour n souches
    # Binariser toutes les variables
    DF_bool = DataFrame(
    age_18_19 = ifelse.(DF.Age .== "18-19",1,0),
    age_20_24 = ifelse.(DF.Age .== "20-24",1,0),
    age_25_34 = ifelse.(DF.Age .== "25-34",1,0),
    age_35_44 = ifelse.(DF.Age .== "35-44",1,0),
    first_time_13_15 = ifelse.(DF.SimulFirstTime .== "13-15",1,0),
    first_time_16_17 = ifelse.(DF.SimulFirstTime .== "16-17",1,0),
    first_time_18_19 = ifelse.(DF.SimulFirstTime .== "18-19",1,0),
    first_time_20 = ifelse.(DF.SimulFirstTime .== "20 et +",1,0),
    condomless_0 = ifelse.(DF.SimulCondom .== "0",1,0),
    condomless_1 = ifelse.(DF.SimulCondom .== "1",1,0),
    condomless_2 = ifelse.(DF.SimulCondom .== "2 et +",1,0),
    partner_0 = ifelse.(DF.SimulPartner .== "0",1,0),
    partner_1 = ifelse.(DF.SimulPartner .== "1",1,0),
    partner_2 = ifelse.(DF.SimulPartner .== "2",1,0),
    partner_3_4 = ifelse.(DF.SimulPartner .== "3-4",1,0),
    partner_5 = ifelse.(DF.SimulPartner .== "5 et +",1,0)  
    )

    # Simulation des vaccins en fonction des covariables
    SimulVaccin = estim_Vaccin_log(DF_bool, Intercept_T)
    DF.SimulVaccin = SimulVaccin
    DF_bool[!, :Vaccin] = SimulVaccin

    # Paramètres de l'âge
    W = estim_coefW(length(souche))

    # Génération des infections par souches
    for j in 1 : length(souche)
        infection = estim_HPV_log(DF_bool, W[j,:], parse(Int64, souche[j]), intercept_infection[j])
        df_HPV[!, "HPV_" * souche[j]] = infection
    end 

    # ----------------------------------------------------------------------------
    # Etape 4 : Fusion des données
    data_simul = hcat(DF,df_HPV)

    # print(first(data_simul,5))

    Y2 = sum(Matrix(data_simul[:,7:10]),dims=2)

    data_simul.Y2 = vec(Y2)

    #Méthode par DEMA
    #jsp si on doit laisser ou non l'intercept
    m = lm(@formula(HPV_18 ~ Age + SimulVaccin + Y2),data_simul) 
    tbl = DataFrame(coeftable(m))
    ic = confint(m)

    # colonne nom et coef
    name_col = names(tbl)[1]
    coef_col = names(tbl)[2]
    idx_Y2  = findfirst(==("Y2"), tbl[!, name_col])
    idx_vac = findfirst(==("SimulVaccin"), tbl[!, name_col])
    push!(biaisDEMA,tbl[idx_Y2, coef_col])
    push!(effet_vaccinDEMA,tbl[idx_vac, coef_col])
    
    #IC DEMA
    vrai_effetVac = -10
    lower = ic[idx_vac, 1]
    upper = ic[idx_vac, 2]
    push!(lowerDEMA, lower)
    push!(upperDEMA, upper)
    if lower ≤ vrai_effetVac ≤ upper
    nbDEMA += 1
    end

    #Méthode par LOLA
    m1 = glm(@formula(HPV_18 ~ Age + SimulVaccin),data_simul, Binomial(),LogLink())
    m2 = glm(@formula(Y2 ~ Age + SimulVaccin),data_simul, Poisson())

    tbl1 = DataFrame(coeftable(m1))
    name_col = names(tbl1)[1]
    coef_col = names(tbl1)[2]
    idx_vac = findfirst(==("SimulVaccin"), tbl1[!, name_col])
    effet_vac_biais = tbl1[idx_vac, coef_col]

    tbl2 = DataFrame(coeftable(m2))
    name_col = names(tbl2)[1]
    coef_col = names(tbl2)[2]
    idx_vac = findfirst(==("SimulVaccin"), tbl2[!, name_col])
    biais = tbl2[idx_vac, coef_col]
    push!(biaisLOLA,biais)
    effet_vac = effet_vac_biais - biais
    push!(effet_vaccinLOLA,effet_vac)

    # IC LOLA
    # T = data_simul.SimulVaccin
    # X = hcat(ones(n), T)
    # mu1_hat = predict(m1)
    # mu2_hat = predict(m2)
    # U1 = (data_simul.HPV_18 .- mu1_hat) .* X
    # U2 = (data_simul.Y2 .- mu2_hat) .* X

    # U = hcat(U1, U2)
    # Ubar = mean(U, dims=1)
    # S = transpose(U .- Ubar) * (U .- Ubar) / n

    # W1 = Diagonal(mu1_hat)
    # W2 = Diagonal(mu2_hat)
    # D1 = transpose(X) * W1 * X / n
    # D2 = transpose(X) * W2 * X / n
    # D = Diagonal([D1[2,2], D2[2,2]])

    # cov_matrix = S[2, 4] / (D1[2,2] * D2[2,2])

    var1 = vcov(m1)[2,2]
    var2 = vcov(m2)[2,2]
    seStar = sqrt(var1 + var2 )#- 2 * cov_matrix)
    lower = effet_vac - 1.96 * seStar
    upper = effet_vac + 1.96 * seStar
    push!(lowerLOLA, lower)
    push!(upperLOLA, upper)
    if lower ≤ vrai_effetVac ≤ upper
    nbLOLA += 1
    end

end

# couvertureDEMA
couvertureDEMA = nbDEMA / M
print(couvertureDEMA)

# couvertureLOLA
couvertureLOLA = nbLOLA / M
print(couvertureLOLA)


# Visualisation

histogram(biaisDEMA,
          bins=50,
          xlabel="Biais DEMA",
          ylabel="Fréquence",
          title="Distribution du biais DEMA",
          color=:lightblue)


histogram(biaisLOLA,
          bins=50,
          xlabel="Biais ETIEVANT",
          ylabel="Fréquence",
          title="Distribution du biais LOLA",
          color=:lightblue)

histogram(effet_vaccinDEMA,
          bins=50,
          xlabel="effetVac DEMA",
          ylabel="Fréquence",
          title="Distribution de l'effet_vac DEMA",
          color=:lightblue)

histogram(effet_vaccinLOLA,
          bins=50,
          xlabel="effet_vac ETIEVANT",
          ylabel="Fréquence",
          title="Distribution de l'effet_vac",
          color=:lightblue)
