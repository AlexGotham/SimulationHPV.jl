# SimulationHPV

Ce package a pour but de simuler un jeu de données représentant des femmes entre 18 et 44 ans qui peuvent être potentiellement affectées par le papillomavirus.

La fonction "datasets(n::Int64, target::Vector{Float64}, souche::Vector{String})" retourne le jeu de données simulé et prend en entrée "n" le nombre d'individus, "target" un vecteur contenant les incidences souhaitées pour chaque infection et "souches" un vecteur contenant le numéro des souches HPV à simuler.

La fonction "borneSomme(vec_inf::Vector{Float64},vec_sup::Vector{Float64})" prend en entrée les bornes inférieures et les bornes supérieures des proportions de la variable que l'on souhaite simuler et retourne un vecteur qui contient les proportions tirées aléatoirement entre les intervalles, avec la condition que la somme des proportions soit égale à 1.

La fonction "estim_Vaccin(df_bool::DataFrame, prop_vaccin::Float64)" prend en entrée un dataframe de valeur binaire représentant les facteurs de confusion et un float contenant la proportion de la population qui est vaccinée. Elle retourne un vecteur binaire indiquant si l'individu est vacciné ou non.

La fontion "sigmoid(x::Float64)" prend en entier un vecteur de float et retourne un vecteur de probabilité.

La fonction "estim_b0(incidence::Float64, reg::Float64)" prend en paramètre l'incidence souhaitée et le résultat de la combinaison linéaire servant à estimer les infections. Et retourner un intercept à ajouter dans la combinaison linéaire afin de contrôler l'incidence.

La fonction "estim_Confusion(mod_var::Vector{String}, prop_var::Vector{Float64}, tauxEvolution::Vector{Vector{Float64}}, n::Int64, proportion_age::Vector{Float64})" prend en entrée le vecteur des modalités de l'habitude sexuelle que l'on souhaite simuler, les proportions par modalités. Les coefficients servant à faire varier les évolutions pour chaque croisement de modalités car on simule les habitudes sexuelles en fonction de l'âge. Le nombre d'individus et les proportions d'individus en fonction de leur âge.

La fonction "estim_HPV(df_bool::DataFrame, W::Vector{Float64}, souche::Int64, target::Float64)" va estimer les infractions aux souches HPV. Elle prend en entrée un dataframe de valeur binaire représentant les facteurs de confusion, les coefficients représentant les effectifs des facteurs de confusion, le nom de la souche et l'incidence.
