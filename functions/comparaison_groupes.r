# Retourne les sous listes de taille n-1 de 'list'
# liste : groupe d'amis (liste d'entiers)
get_sub_lists <- function(liste){
  taille <- length(liste)
  res <- vector(mode= "list", length= taille-1)
  
  for(i in 1:taille){
    sous_liste <- vector(mode= "list", length= taille-1)
    for(j in 1:taille){
      if(i != j){
        sous_liste <- append(sous_liste, liste[[j]])
      }
    }
    res <- append(res, sous_liste)
    print(unlist(sous_liste))
  }
  
  return (res)
}


# Retourne TRUE si la liste est contenue dans la liste de listes
is_list_of <- function(list, listOfLists){
  list_to_compare <- sort(unlist(list))
  for(list_test in listOfLists){
    list_to_test <- sort(unlist(list_test))
    if(identical(list_to_compare, list_to_test)) return (TRUE)
  }
  return (FALSE)
}

# Tri la liste d'amis par taille de groupe decroissante
# amis : ensemble des groupes d'amis
tri_amis_decroissant <- function(amis){
  res <- list()
  taille <- length(amis)
  tailles_groupes <- c() #liste des tailles de groupes d'amis
  for(i in 1:taille){
    tailles_groupes <- append(tailles_groupes, length(amis[[i]]))
  }
  
  #tri des groupes d'amis selon leurs tailles
  indices_tri <- sort(tailles_groupes, decreasing= TRUE, index.return= 1)$ix #ix donne les indices

  res <- amis[indices_tri]
  
  return (res)
}

#   Idee :
#   Parcourir les groupes de maniere decroissante
#   Pour chaque groupe, regarder s'il est present dans groupe 2 tel quel
#   Si c'est le cas, c'est le plus grand groupe
#   Sinon on teste avec un sous-echantillon de taille n-1
#   On compare ensuite les tailles des meilleurs echantillons retenus pour chaque groupe

# Recupere les proteines formant le plus grand groupe se retrouvant toujours ensemble parmi amis1 et amis2
# amis1 : ensemble des groupes d'amis (liste de listes d'entiers)
# amis2 : ensemble des groupes d'amis (liste de listes d'entiers)
get_solid_elements <- function(amis1, amis2){
  taille1 <- length(amis1)
  taille2 <- length(amis2)
  
  amis1 <- tri_amis_decroissant(amis1)
  amis2 <- tri_amis_decroissant(amis2)
  
  for(ami1 in 1:taille1) #parcours des groupes d'amis de amis1
  {
    if(is_list_of(amis1[[ami1]], amis2)){
      print("groupe equivalent trouve !")
    }
  }

  
    
}





