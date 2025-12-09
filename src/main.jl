# Fonction générative du jeu de données
using DataFrames
using Distributions

include("main.jl")
include("estim_b0.jl")
include("estim_Confusion.jl")
include("estim_HPV.jl")
include("estim_proportionIC.jl")
include("estim_Sigmoid.jl")
include("estim_Vaccin.jl")

function datasets(n::Int64, target::Vector{Float64}, souche::Vector{String})
    # Initialisation 
    # ----------------------------------------------------------------------------
    # Etape 1 : Simulation des variables
    # simulation age
    modalite_age = ["18-19","20-24","25-34","35-44"]
    born_inf_age = [5,16,34.7,37.1]/100
    born_sup_age = [6.8,19.2,39.2, 42.3]/100
    prop_age = borneSomme(born_inf_age,born_sup_age)
    age = wsample(modalite_age, prop_age, n) 
    
    # simulation  HPV vaccination status
    modalite_vac = ["Vaccinated","Unvaccinated"]
    prop_vac = rand(Uniform(50,68))/100
    
    # simulation age first sex
    modalite_age_sex = ["13-15","16-17","18-19","20 et +"]
    born_inf_age_sex = [24.1,40.7,17.3,9.6]/100
    born_sup_age_sex = [28.1,45.8,21.6,13.3]/100
    prop_age_sex = borneSomme(born_inf_age_sex,born_sup_age_sex)
    
    # proportion condomless sex with partner, past year
    modalite_condomless = ["0","1","2 et +"]
    born_inf_condomless = [16.5,70.9,7.2]/100
    born_sup_condomless = [20.7,75.5,9.4]/100
    prop_condomless = borneSomme(born_inf_condomless,born_sup_condomless)

    # proportion de Total partners, past 5 years
    modalite_partner = ["0","1","2","3-4","5 et +"]
    born_inf_partner = [0.9,54.9,12.7,11.4,12.8]/100
    born_sup_partner = [2.4,59.8,15.8,14.4,15.6]/100
    prop_partner = borneSomme(born_inf_partner,born_sup_partner)
    # ----------------------------------------------------------------------------
    # Etape 2 : Simulation des variables avec interaction de l'âge
    # Première fois en fonction de l'âge
    taux_age_sex = [
        [1.2, 1.6, 1.2, 0], 
        [1,1.6,1,0.4],  
        [1,1.4,1,0.6],   
        [1,1.4,1,0.6]
    ]
    SimulFirstTime = estim_Confusion(modalite_age, prop_age, taux_age_sex , n, prop_age)

    # Utilisation de protection en fonction de l'âge
    taux_condomless = [
        [1.2, 1.6, 1.2], 
        [1,1.6,1],  
        [1,1.4,1],   
        [1,1.4,1]
    ]
    SimulCondom = estim_Confusion(modalite_condomless, prop_condomless, taux_condomless, n, prop_age)

    # Nombre de partenaire en fonction de l'âge
    taux_partner = [
        [0, 2, 1.2,1,0.8], 
        [0.3,1.8,1,1,0.9],  
        [0.3,1.5,1.1,1.1,1],   
        [0.3,1.7,1.1,1.1,0.8]
    ]
    SimulPartner = estim_Confusion(modalite_partner, prop_partner, taux_partner, n, prop_age)

    # Regroupement des données dans un dataframe
    Age = sort(age)
    DF = DataFrame(Age = Age,  SimulFirstTime = SimulFirstTime, SimulCondom = SimulCondom, SimulPartner = SimulPartner)

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
    SimulVaccin = estim_Vaccin(DF_bool, prop_vac)
    DF_bool[!, :Vaccin] = SimulVaccin

    # Paramètres de l'âge
    W = [
        [0.5, 1, 1, -0.5], 
        [1.5,1,-0.5,1],  
        [2.5,-0.5,1,2.5],  
        [0.5,1,0.5,2],   
        [2,1.5,-0.5,2.5]
    ]

    # Génération des infections par souches
    df_infection = DataFrame()
    for i in 1 : length(souche)
        infection = estim_HPV(DF_bool, W[i], parse(Int64, souche[i]), target[i])
        df_infection[!, "HPV_" * souche[i]] = infection
    end 

    # ----------------------------------------------------------------------------
    # Etape 4 : Fusion des données
    DATA = hcat(DF,df_infection)
    return(DATA)
end

