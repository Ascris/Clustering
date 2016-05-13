#retourne vrai si les proteines presentes aux indices "triplet" sont triangulaires
is_triangular <- function(matDist, triplet){
  #regarder si D_p1p2 <= Dp1_p3 ... etc
  i <- triplet[1]
  j <- triplet[2]
  k <- triplet[3]
  
  return((matDist[i,j] <= matDist[i,k]+matDist[j,k]) &&
     (matDist[i,k] <= matDist[i,j]+matDist[j,k]) &&
     (matDist[j,k] <= matDist[i,j]+matDist[i,k]))
}

#calcul de la triangulation de 3 points au hasard parmi ceux de matDist
#matDist = matrice de distances etudiee
#nb_prot = nombre de proteines du fichier
#nb_ite = nombre d'iterations souhaite
#retourne le pourcentage de "triplets gagnants"
triangular <- function(matDist, nb_prot, nb_ite){
  nb_triangulation <- 0
  for(z in 1:nb_ite){
    triplet <- sample(1:nb_prot, 3)
    triangulation_OK <- is_triangular(matDist, triplet)
    if(triangulation_OK) nb_triangulation <- nb_triangulation + 1
  }
  triplets_gagnants <- nb_triangulation/nb_ite*100
  return (triplets_gagnants)
}

#retourne vrai si les proteines presentes aux indices "triplet" sont triangulaires
is_triangular_ultrametric <- function(matDist, triplet){
  #regarder si D_p1p2 < Dp1_p3 ... etc
  i <- triplet[1]
  j <- triplet[2]
  k <- triplet[3]
  
#   print(paste("matDist[i,j] = ", matDist[i,j]))
#   print(paste("matDist[i,k] = ", matDist[i,k]))
#   print(paste("matDist[j,k] = ", matDist[j,k]))
#   
#   print(paste("res = ", ((matDist[i,j] <= max(matDist[i,k], matDist[j,k])) &&
#                            (matDist[i,k] <= max(matDist[i,j], matDist[j,k])) &&
#                            (matDist[j,k] <= max(matDist[i,j]+matDist[i,k])))))
  
  return((matDist[i,j] <= max(matDist[i,k], matDist[j,k])) &&
           (matDist[i,k] <= max(matDist[i,j], matDist[j,k])) &&
           (matDist[j,k] <= max(matDist[i,j], matDist[i,k])))
}

#calcul de la triangulation de 3 points au hasard parmi ceux de matDist
#matDist = matrice de distances etudiee
#nb_prot = nombre de proteines du fichier
#nb_ite = nombre d'iterations souhaite
#retourne le pourcentage de "triplets gagnants"
triangular_ultrametric <- function(matDist, nb_prot, nb_ite){
  nb_triangulation <- 0
  for(z in 1:nb_ite){
    triplet <- sample(1:nb_prot, 3)
    triangulation_OK <- is_triangular_ultrametric(matDist, triplet)
    if(triangulation_OK) nb_triangulation <- nb_triangulation + 1
  }
  triplets_gagnants <- nb_triangulation/nb_ite*100
  return (triplets_gagnants)
}