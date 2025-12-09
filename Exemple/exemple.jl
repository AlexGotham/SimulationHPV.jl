using SimulationHPV

# Fonction sigmoïde
print(sigmoid.([12.0,12.0,130]))

# Fonction Confusion
taux_age_sex = [[0.01,0.03 , 0.02, 0],[1,1.6,1,0.4],[1,1.4,1,0.6],[1,1.4,1,0.6]]
propo_age = [0.06,0.18,0.37,0.39]
moda_ageSexe = ["18-19","20-24","25-34","35-44"]
propo_ageSexe = [0.1,0.2,0.4,0.3]
n = 5
estim_Confusion(moda_ageSexe,propo_ageSexe,taux_age_sex,n,propo_age)

# Fonction Estim_Vaccin
df_bool_vaccin = DataFrame(rand(Bool, 100, 12),:auto)
prop_vaccin_simul = 0.55
print(estim_Vaccin(df_bool, prop_vaccin_simul))

# Fonction Estim_b0
print(estim_b0.([0.05,0.06],[12.9,13]))

# Fonction Estim_HPV
prop = borneSomme([16.5,70.9,7.2]/100,[20.7,75.5,9.4]/100)
df_bool_HPV = DataFrame(rand(Bool, 100, 15),:auto)
print(estim_HPV(df_bool_HPV,[0.4,0.2],6,0.05))

# Fonction datasets
n = 1000
target = [0.22,0.055,0.03,0.05,0.025]
souche = ["6","11","16","18","31"]
data_simul = datasets(n,target,souche)
print(first(data_simul,5))

