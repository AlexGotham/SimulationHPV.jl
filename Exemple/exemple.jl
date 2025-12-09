using SimulationHPV

n = 1000
target = [0.22,0.055,0.03,0.05,0.025]
souche = ["6","11","16","18","31"]
data_simul = datasets(n,target,souche)
print(first(data_simul,5))