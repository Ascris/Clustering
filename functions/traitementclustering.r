list.of.packages <- c("cluster")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(cluster)

#' Teste le clustering par 2 à nb_clust (>2) et retourne le cas où la moyenne etait maximal
#'
#' @param mat : matrice de distance ou de robustesse a tester (matrice de reels ou d'entiers)
#' @param nb_clust : nombre de groupes maximum a former avec PAM (entier)
#'
#' @return meilleur nombre de groupes a faire, selon la moyenne des silhouettes (entier)
#' @export Teste le clustering par 2 à nb_clust (>2) et retourne le cas où la moyenne etait maximal
#'
#' @examples nb_clusters <- get_best_clustering(data_1751, 5)
get_best_clustering <- function(mat, nb_clust){
  if(nb_clust <= 2) return (2)
  res <- 2
  val_max <- pam(mat, 2)$silinfo$avg.width
  for(i in 3:nb_clust){
    tmp_max <- pam(mat, i)$silinfo$avg.width
    if(tmp_max > val_max){
      val_max <- tmp_max
      res <- i
    }
  }
  return (res)
}

#' Compte le nombre de groupes differents presents dans 'groupes'
#'
#' @param groupes : vecteur d'entiers (resultat de pam$clustering par exemple)
#'
#' @return nombre de groupes differents (entier)
#' @export Compte le nombre de groupes differents presents dans 'groupes'
#'
#' @examples diff_groups <- getNbGroups(myPam$clustering)
getNbGroups <- function(groupes){
  nb_groups <- 0
  checked_groups <- c()
  for(i in 1:length(groupes)){
    if(groupes[i] %in% checked_groups){
      #ne rien faire
    } else {
      checked_groups <- append(checked_groups, groupes[i])
      nb_groups <- nb_groups + 1
    }
  }
  return (nb_groups)
}

###################################
######## ECHANTILLONNAGE ##########
###################################

#' A partir de 'data1751', releve 'taille_ech' individus, y applique PAM puis stocke les groupes dans des fichiers
#'
#' @param data_X : matrice de distance des X proteines (matrice de reels)
#' @param names_prot : noms des proteines (vecteur de chaines de caracteres)
#' @param pourcentage : pourcentage de l'echantillon a prelever (entier)
#' @param Cl : nombre de clusters a former dans le PAM (entier)
#'
#' @return moyenne des silhouettes retournees par PAM
#' @export A partir de 'data1751', releve 'taille_ech' individus, y applique PAM puis stocke les groupes dans des fichiers
#'
#' @examples res <- echantillonnage(data1751, names1751, 80, 5)
echantillonnage <- function(data_X, names_prot, pourcentage, Cl){
  nombre_proteines <- length(data_X[1,])
  
  rd <- sort(sample(1:nombre_proteines, nombre_proteines*pourcentage/100)) #proteines selectionnees
  mat <- data_1751[rownames=rd, colnames=rd] #extraction sous-matrice
  
  pamEchantillon <- pam(mat, Cl)
  plot(pamEchantillon)
  
  groupes <- pamEchantillon$clustering #recuperation clusters
  ecrire_fic_groupe(groupes, names_prot, rd) #stockage des clusters dans des fichiers
  
  return (pamEchantillon$silinfo$avg.width)
}

#' Effectue la moyenne des resultats de N echantillonnages de X individus de data_X
#'
#' @param data_X : matrice de distance de X proteines (matrice de reels)
#' @param names_prot : vecteur des noms des proteines (vecteur de chaines de caracteres)
#' @param X : taille de l'echantillon a prelever (entier)
#' @param Cl : nombre de clusters a former dans le PAM (entier)
#' @param N : nombre d'echantillonnages a effectuer (entier)
#'
#' @return moyenne des echantillonnages (entier)
#' @export Effectue la moyenne des resultats de N echantillonnages de X individus de data_X
#'
#' @examples moy <- moyenne_echantillonnage(data_403, names403, 80, nb_clusters, 10)
moyenne_echantillonnage <- function(data_X, names_prot, X, Cl, N){
  tot <- 0.0
  for(i in 1:N){
    tot <- tot + echantillonnage(data_X, names_prot, X, Cl)
  }
  moy_ech <- tot/N
  return(moy_ech)
}

###############################################
############ COMPTAGE CRITERE PAM #############
###############################################

#' Dans le cluster 'num_groupe', compte combien d'elements respectent le critere
#'
#' @param pam : resultat du PAM qui nous interesse (objet de classe 'pam')
#' @param namesX  : noms des X proteines placees dans les clusters (vecteur de chaines de caracteres)
#' @param num_groupe : numero du groupe interesse (entier)
#' @param critere : comptage des proteines selon ce critere (Thermophile, Bacteria,..etc) (chaine de caracteres)
#'
#' @return nombre d'elements respectant le critere (entier)
#' @export Dans le cluster 'num_groupe', compte combien d'elements respectent le critere
#'
#' @examples count_in_groupe(pam$clustering, 1, "Bacteria")
count_in_group <- function(pam, namesX, num_groupe, critere){
  groupes <- pam$clustering #liste des individus repartis dans les clusters
  nb_individus <- pam$clusinfo[num_groupe, 1] #taille groupe 'num_groupe'
  res <- 0 #proteines qui respectent le critere
  type_critere <- which_critere(critere) #critere = regne(A,B ou E) | type(TOP ou RG) | milieu(T, M, P ou H)
  
  for(i in 1:length(groupes)){
    if(num_groupe == groupes[i]) #cluster qui nous interesse
    {
      corresponding <- switch(type_critere,
                              "regne"= testRegne(namesX[i], critere),
                              "type"= testType(namesX[i], critere),
                              "milieu"= testMilieu(namesX[i], critere))
      res <- res + corresponding
    }
  }
  print(paste("groupe ", num_groupe, " - ", res, " / ", nb_individus))
  return (res)
}

#' Retourne le nombre total d'individus respectant le critere
#'
#' @param pam : resultat de la fonction 'pam' (objet de la classe 'pam')
#' @param namesX : noms des proteines presentes dans le PAM (vecteur de chaines de caracteres)
#' @param critere : "RG", "A"...etc (chaine de caracteres)
#'
#' @return nombre total d'individus respectant le critere (entier)
#' @export Retourne le nombre total d'individus respectant le critere
#'
#' @examples nb_ind <- count_in_all_groups(pam1751, names1751, "M")
count_in_all_groups <- function(pam, namesX, critere){
  nb_individus <- length(pam$clustering)
  res <- 0
  diff_groups <- getNbGroups(pam$clustering)
  for(i in 1:diff_groups){
    res <- res + count_in_group(pam, namesX, i, critere)
  }
  print(paste("Total - ", get_full_critere(critere), " : ", res, " / ", nb_individus))
  return (res)
}

#' Verifie que la proteine respecte bien le critere
#'
#' @param proteine : indice de la proteine (entier)
#' @param critere : "TOP", "A", ...etc (chaine de caracteres)
#'
#' @return TRUE si la proteine respecte le critere
#' @export Verifie que la proteine respecte bien le critere
#'
#' @examples valid <- est_valide("AcamarTOPa_M_B", "T")
est_valide <- function(proteine, critere){
  type_critere <- which_critere(critere) #critere = regne(A,B ou E) | type(TOP ou RG) | milieu(T, M, P ou H)
  corresponding <- switch(type_critere,
                          "regne"=testType(proteine, critere),
                          "type"=testParticularite(proteine, critere),
                          "milieu"=testMilieu(proteine, critere))
  if(1 == corresponding) return (TRUE)
  else return (FALSE)
}

#' A partir du pam et des noms des proteines presentes dedans, compte combien respectent les criteres
#'
#' @param pam : resultat de la fonction 'pam' (objet de classe 'pam')
#' @param namesX : noms des proteines presentes dans le pam (vecteur de chaines de caracteres)
#' @param criteres : vecteur de criteres de type "RG", "A" ..etc (chaine de caracteres)
#'
#' @return nombre de proteines respectant tous les criteres
#' @export A partir du pam et des noms des proteines presentes dedans, compte combien respectent les criteres
#'
#' @examples nb_ind <- count_multi_criteres(pam1751, names1751, c("M","RG"))
count_multi_criteres <- function(pam, namesX, criteres){
  groupes <- pam$clustering #liste des individus repartis dans les clusters
  res <- 0 #nombre de proteines qui respectent le critere
  
  for(i in 1:length(groupes)){
    nom_prot <- namesX[i] #proteine courante
    all_criteres_valides <- 0 #nombre de criteres valides
    for(j in 1:length(criteres)){
      if(est_valide(nom_prot, criteres[j])){
        all_criteres_valides <- all_criteres_valides + 1
      }
    }
   if(all_criteres_valides == length(criteres)) res <- res + 1 #tous les criteres sont valides
  }
  return (res)
}


###############################################
########## FONCTIONS SUR LES GROUPES ##########
###############################################

#' Verifie si la liste est contenue dans la liste de listes
#'
#' @param list : liste a tester
#' @param listOfLists : liste de listes
#'
#' @return TRUE si la liste est contenue dans la liste de listes (booleen)
#' @export Verifie si la liste est contenue dans la liste de listes
#'
#' @examples est_contenue_dans <- is_list_of(list(1,2), list(list(1,2), list(4,5)))
is_list_of <- function(list, listOfLists){
  list_to_compare <- sort(unlist(list))
  for(list_test in listOfLists){
    list_to_test <- sort(unlist(list_test))
    if(identical(list_to_compare, list_to_test)) return (TRUE)
  }
  return (FALSE)
}

#' Retourne la liste des amis 'forts' et DIRECTS de x dans la matrice de robustesse
#'
#' @param matRobustesse : matrice de robustesse (matrice d'entiers)
#' @param x : indice de la proteine dont on veut les amis (entier)
#' @param force : puissance du lien qui doit lier les proteines au sein d'un groupe (entier)
#'
#' @return liste des amis 'forts' et DIRECTS de x
#' @export Retourne la liste des amis 'forts' et DIRECTS de x dans la matrice de robustesse
#'
#' @examples amis_de_x <- get_friends(matRobustesse, 2, 19)
get_friends <- function(matRobustesse, x, force){
  friends <- c()
  taille <- length(matRobustesse[1,])
  for(i in 1:taille)
  {
    if(force <= matRobustesse[x,i]) friends <- union(friends, i)
  }
  return (friends)
}

#' Recupere le reseau d'amis
#'
#' @param matRobustesse : matrice de robustesse (matrice d'entiers)
#' @param res : groupe en construction (liste d'entiers)
#' @param liste_amis : amis voulant entrer dans le groupe en formation (liste d'entiers)
#' @param force : puissance du lien qui doit lier les proteines au sein d'un groupe (entier)
#'
#' @return reseau d'amis des individus de 'liste_amis'
#' @export Recupere le reseau d'amis
#'
#' @examples reseau_ami <- get_res(matRobustesse, local_list, friends_i, force)
get_res <- function(matRobustesse, res, liste_amis, force){
  if(is_list_of(liste_amis, res)) return (res)
  
  for(ami in liste_amis){
    if(ami %in% res){
      #rien ne se passe
    } else {
        res <- union(res, ami) #ajout de l'ami courant a la liste s'il n'y est pas deja
        res <- union(res, get_res(matRobustesse, res, get_friends(matRobustesse, ami, force), force))
    }
  }
  return (res)
}

#' Recupere la liste des reseaux d'amis, 2 proteines sont amies lorsqu'elle ont un lien fort (de valeur >= 'force')
#'
#' @param matRobustesse : matrice de robustesse (matrice d'entiers)
#' @param force : puissance du lien qui doit lier les proteines au sein d'un groupe (entier)
#'
#' @return liste des reseaux d'amis
#' @export Recupere la liste des reseaux d'amis, 2 proteines sont amies lorsqu'elle ont un lien fort (de valeur >= 'force')
#'
#' @examples liste_de_reseaux <- get_all_friends(matRobustesse, 21)
get_all_friends <- function(matRobustesse, force){
  liste_liste_amis <- list() #ensemble des reseaux d'amis
  taille <- length(matRobustesse[1,])
  force_max <- matRobustesse[1,1]
  for(i in 1:taille) #i = proteine courante
  {
    print(paste("get friends of protein ", i, "/", taille, sep=""))
    local_list <- list()
    friends_i <- get_friends(matRobustesse, i, force) #amies de proteine i
    if(force == force_max) #le lien est maximal
    {
      friends_i <- get_friends(matRobustesse, i, force) #amies de proteine i
      liste_liste_amis[[i]] <- friends_i
    } else #le lien n'est pas maximal
    {
      local_list <- get_res(matRobustesse, local_list, friends_i, force) #reseau d'amies de la proteine i
      liste_liste_amis[[i]] <- local_list
    }
  }
  
  # enleve les doublons de liste_liste_amis et les groupes de taille 1 (singletons)
  real_list <- list()
  real_list <- suppression_doublons(liste_liste_amis)
  
  return (real_list)
}

#' Retourne les id des proteines de 1 a nb_prot qui ne sont dans aucun reseau d'amis de amisX
#'
#' @param amisX : ensemble des reseaux d'amis (liste de liste d'entiers)
#' @param nb_prot : nombre total de proteines dans le jeu de donnees (entier)
#'
#' @return id des proteines de 1 a nb_prot qui ne sont dans aucun reseau d'amis de amisX
#' @export Retourne les id des proteines de 1 a nb_prot qui ne sont dans aucun reseau d'amis de amisX
#'
#' @examples singletons <- find_singletons(amis1751, 1751)
find_singletons <- function(amisX, nb_prot){
  all_friends <- c() # liste des proteines presentes dans les réseaux d'amis
  for(i in 1:length(amisX)){
    for(j in 1:length(amisX[[i]])){
      all_friends <- append(all_friends, amisX[[i]][[j]])
    }
  }
  sort(all_friends)
  singletons <- setdiff(1:nb_prot, all_friends)
  return (singletons)
}

###############################################
######### FONCTIONS MATRICE ROBUSTESSE ########
###############################################

#' Applique PAM a 'dataX' avec des clusters de 2 a 'max_clusters' puis stocke les resultats
#' dans une matrice[X,X] dans laquelle la case [i,j] correspond a la robustesse du lien entre
#' i et j (nombre de fois que i et j se sont retrouves dans le meme cluster)
#'
#' @param dataX : matrice de distance (matrice de reels)
#' @param taille : dimension de la matrice de distance (entier)
#' @param min_clusters : nombre minimal de clusters a tester (entier)
#' @param max_clusters : nombre maximal de clusters a tester (entier)
#'
#' @return matrice de robustesse
#' @export Applique max-min fois PAM et construit la matrice de robustesse
#'
#' @examples matRob <- build_mat_rob(data403, 403, 5, 25)
build_mat_rob <- function(dataX, taille, min_clusters, max_clusters){
  matRob <- matrix(nrow= taille, ncol=taille, data = rep(0,taille*taille))
  for(c in min_clusters:max_clusters){
    current_pam <- pam(dataX, c)
    groupes <- current_pam$clustering #vecteur de 'taille' entiers
    print(paste("construction pam : ", c, " sur ", max_clusters))
    for(i in 1:taille){
      for(j in 1:i){
        if(groupes[i] == groupes[j])
        {
          matRob[i,j] <- matRob[i,j] + 1
          if(i != j) matRob[j,i] <- matRob[j,i] + 1
        }
      }
    }
  }
  return (matRob)
}


#' Construit une matrice de robustesse basee sur dataX de taille m.
#' Cette derniere contient le nombre de fois que 2 proteines ont ete classees ensemble par PAM
#'
#' @param dataX : matrice de distance (matrice de reels)
#' @param m : taille de la matrice de robustesse a creer (entier)
#' @param k : nombre de groupes demandes a PAM (entier)
#' @param n : nombre d'iterations de l'algo PAM a effectuer (entier)
#'
#' @return matrice de robustesse de taille m
#' @export Construit une matrice de robustesse basee sur dataX de taille m
#'
#' @examples matRob <- build_k_mat_rob(data403, 403, 5, 1)
build_k_mat_rob <- function(dataX, m, k, n){
  matRob <- matrix(nrow= m, ncol=m, data = rep(0,m)) #initilisation a 0 partout
  for(taille in 1:n){
    current_pam <- pam(dataX, k)
    groupes <- current_pam$clustering
    for(i in 1:m){
      for(j in 1:i){
        if(groupes[i] == groupes[j]) #les proteines i et j ont ete classees ensemble
        {
          matRob[i,j] <- matRob[i,j] + 1
          if(i != j) matRob[j,i] <- matRob[j,i] + 1 #remplissage de la diagonale superieure
        }
      }
    }
  }
  return (matRob)
}


#' Construit une matrice de robustesse basee sur dataX de taille m en testant un PAM de min a max groupes
#'
#' @param dataX : matrice de distance (matrice de reels)
#' @param m : taille de la matrice de robustesse a creer (entier)
#' @param n : nombre d'iterations de l'algo PAM (entier)
#' @param min : minimum de groupes voulus avec PAM (entier)
#' @param max : maximum de groupes voulus avec PAM (entier)
#'
#' @return matrice de robustesse basee sur dataX de taille m en testant un PAM de min a max groupes
#' @export Construit une matrice de robustesse basee sur dataX de taille m en testant un PAM de min a max groupes
#'
#' @examples matRob <- build_all_mat_rob(data403, 403, 1, 2, 7)
build_all_mat_rob <- function(dataX, m, n, min, max){
  super_matRob <- matrix(nrow= m, ncol=m, data = rep(0,m)) #initialisation a 0
  
  for(k in min:max){
    print(paste("Construction matRob ", k))
    current_matRob <- build_k_mat_rob(dataX, m, k, n)
    
    #sauvegarde les matRob dans le dossier matRob
    fic_name <- paste("matRob", m, "_", k, sep="")
    save_matrix("~/R/resultats/matRob/", fic_name, current_matRob)
    
    super_matRob <- super_matRob + current_matRob
  }
  
  #sauvegarde la super_matRob dans le dossier matRob
  fic_name <- paste("SUPER_matRob", m, "_", k, sep="")
  save_matrix("~/R/resultats/matRob/", fic_name, super_matRob)
  
  return (super_matRob)
}


#' A partir des donnees, du nombre min et max de clusters a tester ainsi que des noms des proteines,
#' cree la super matrice de robustesse, les cliques regne-type et regne-milieu et le fichier contenant
#' les groupes de proteines ayant un lien fort entre elles
#'
#' @param data_X : matrice de distance (matrice de reels)
#' @param X : nombre de proteines presentes (entier)
#' @param n : nombre de fois que PAM a ete repete pour definir les groupes (entier)
#' @param min : nombre minimum de clusters a former (entier)
#' @param max : nombre maximum de clusters a former (entier)
#' @param namesX : noms des proteines (vecteur de chaines de caracteres)
#' @param dim : valeur de diminution jusqu'a laquelle on va tester la robustesse des groupes (entier)
#'
#' @return la super matrice de robustesse
#' @export Cree la super matrice de robustesse et les cliques associees
#'
#' @examples stabilite_groupes403 <- stabilite_groupes(data_403, 403, 1, 5, 25, names403, 5)
stabilite_groupes <- function(data_X, X, n, min, max, namesX, dim){
  super_matRob <- build_all_mat_rob(data_X, X, n, min, max)
  for(i in 0:dim){
    robustesse <- (n*length(min:max))-i #seuil de robustesse devant lier les proteines entre elles
    
    criteres <- c("TOP", "T")
    nom_clique_1 <- paste("clique", X, "_", criteres[1], "_n", robustesse, ".png", sep="")
    nom_clique_2 <- paste("clique", X, "_", criteres[2], "_n", robustesse, ".png", sep="")
    
    
    nom_fic_sommets <- paste("proteines", X, ".txt", sep="")
    nom_fic_arcs <- paste("robustesse", X, ".txt", sep="")
    creer_fichier_sommets(nom_fic_sommets, namesX)
    creer_fichier_arcs(nom_fic_arcs, super_matRob, namesX, robustesse)
    
    super_clique_1 <- create_clique(nom_clique_1, nom_fic_sommets, nom_fic_arcs, criteres[1], X)
    super_clique_2 <- create_clique(nom_clique_2, nom_fic_sommets, nom_fic_arcs, criteres[2], X)
    
    super_amis <- get_all_friends(super_matRob, robustesse) #recherche du reseau d'amis
    
    nom_fic_groupes <- paste("super_amis", X, "_n", robustesse, ".txt", sep="")
    ecriture_fichier_groupes("clique", nom_fic_groupes, super_amis, namesX)
  }
   return (super_matRob)
}

#' Supprime les groupes d'amis deja presents dans la liste d'amis (liste de liste d'entiers)
#'
#' @param liste_amis : liste dans laquelle enlever les doublons
#'
#' @return liste de liste d'amis sans les doublons
#' @export Supprime les groupes d'amis deja presents dans la liste d'amis (liste de liste d'entiers)
#'
#' @examples real_list <- suppression_doublons(amis403)
suppression_doublons <- function(liste_amis){
  final_list <- list()
  for(i in 1:length(liste_amis)){
    if(is_list_of(liste_amis[i], final_list)){
      # print("liste deja presente !")
    } else {
      if(1 >= length(liste_amis[[i]])){
        # print("singleton non accepte !")
      } else {
        # print("ajout OK ")
        final_list[length(final_list)+1] <- liste_amis[i]
      }
    }
  }
  return (final_list)
}

#' Retournement de la matrice selon le X donne
#'
#' @param matRob : matrice de robustesse a inverser (matrice d'entiers)
#' @param X : valeur d'inversement de la matrice (entier)
#'
#' @return matrice de robustesse inversee
#' @export Retournement de la matrice selon le X donne
#'
#' @examples matRob <- retournementMat(matRob, 19)
retournementMat <- function(matRob, X){
  for(i in 1:length(matRob[1,])){
    for(j in 1:length(matRob[1,])){
      matRob[i, j] <- X - matRob[i, j]
    }
  }
  return (matRob)
}

#' A partir de la coupe faite sur le resultat d'un hclust, retourne une liste de liste d'amis
#'
#' @param coupe : resultat de la fonction cutree appliquee au resultat d'un hclust (vecteur d'entiers)
#' @param taille : groupes classes ensemble par le cutree (liste de listes d'amis)
#'
#' @return liste de liste d'amis associee a la coupe par hclust et cutree
#' @export A partir de la coupe faite sur le resultat d'un hclust, retourne une liste de liste d'amis
#'
#' @examples amis_coupe10 <- build_friend_list(coupe_10, 10)
build_friend_list <- function(coupe, taille){
  taille_coupe <- length(coupe)
  friend_list <- vector(mode= "list", length= taille)
  for(i in 1:taille_coupe){
    groupe <- coupe[i]
    friend_list[[groupe]] <- union(friend_list[[groupe]], i)
  }
  
  return (friend_list)
}

#' Trie la liste d'amis sur leurs valeurs decroissantes
#'
#' @param amis : liste de listes d'amis
#'
#' @return liste de listes d'amis triee
#' @export Trie la liste d'amis sur leurs valeurs decroissantes
#'
#' @examples amis403 <- tri_amis_decroissant(amis403)
tri_amis_decroissant <- function(amis){
  res <- list()
  taille <- length(amis)
  tailles_groupes <- c() #liste des tailles de groupes d'amis
  for(i in 1:taille){
    tailles_groupes <- append(tailles_groupes, length(amis[[i]]))
  }
  
  #tri des groupes d'amis selon leurs tailles
  indices_tri <- sort(tailles_groupes, decreasing= TRUE, index.return= 1)$ix #ix donne les indices
  
  res <- amis[indices_tri] #tri des amis selon leur taille
  
  return (res)
}

#' A partir des hauteurs de tous les noeuds, retourne le vecteur des distances inter-noeuds
#'
#' @param hauteur_noeuds : hauteurs des noeuds du hclust vecteur de reels
#'
#' @return vecteur des distances inter-noeuds
#' @export A partir des hauteurs de tous les noeuds, retourne le vecteur des distances inter-noeuds
#'
#' @examples distInterNoeuds <- getDistInterNoeuds(hauteur_noeuds)
getDistInterNoeuds <- function(hauteur_noeuds){
  ordre <- sort(hauteur_noeuds, decreasing= TRUE, index.return= 1)$ix
  hauteur_noeuds <- sort(hauteur_noeuds, decreasing= TRUE) #tri necessaire !
  res <- vector(length= length(hauteur_noeuds))
  for(i in 1:(length(hauteur_noeuds)-1)) #pour tous les noeuds qui ont un noeud suivant
  {
    res[i] <- abs(hauteur_noeuds[i] - hauteur_noeuds[i+1])
  }
  res[length(hauteur_noeuds)] <- 0 #dernier noeud sans noeud suivant
  return (res[ordre])
}

#' Retourne le vecteur des hauteurs (distances de la racine) de chaque inter-noeuds
#'
#' @param hauteur_noeuds : hauteurs des noeuds du hclust (vecteur de reels)
#'
#' @return vecteur des hauteurs (distances de la racine) de chaque inter-noeuds
#' @export Retourne le vecteur des hauteurs (distances de la racine) de chaque inter-noeuds
#'
#' @examples distRacineInterNoeuds <- getDistRacineInterNoeuds(hauteur_noeuds)
getDistRacineInterNoeuds <- function(hauteur_noeuds){
  ordre_initial <- sort(hauteur_noeuds, decreasing= TRUE, index.return= 1)$ix
  hauteur_noeuds <- sort(hauteur_noeuds, decreasing= TRUE) #tri necessaire !
  res <- vector(length= length(hauteur_noeuds))
  for(i in 1:(length(hauteur_noeuds)-1)) #pour tous les noeuds qui ont un noeud suivant
  {
    res[i] <- min(hauteur_noeuds[i], hauteur_noeuds[i+1]) + abs((hauteur_noeuds[i] - hauteur_noeuds[i+1])/2)
  }
  res[length(hauteur_noeuds)] <- 0 #dernier noeud sans noeud suivant
  return (res[ordre_initial])
}


#' Retourne le classement des individus maximisant la distance inter-noeuds et la hauteur
#'
#' @param distIN : distance inter-noeuds (vecteur de reels)
#' @param distRacineIN : hauteur des milieux d'inter-noeuds
#'
#' @return classement des individus maximisant la distance inter-noeuds et la hauteur
#' @export Retourne le classement des individus maximisant la distance inter-noeuds et la hauteur
#'
#' @examples classementNoeuds <- getBestNodes(distInterNoeuds, distRacineInterNoeuds)
getBestNodes <- function(distIN, distRacineIN){
  res <- c(rep(0, length(distIN)))
  
  ordre_distIN <- sort(distIN, decreasing= TRUE, index.return= 1)$ix
  ordre_distRacineIN <- sort(distRacineIN, decreasing= TRUE, index.return= 1)$ix
  
  #classement des distances inter-noeuds
  for(i in 1:length(ordre_distIN)){
    res[ordre_distIN[i]] <- res[ordre_distIN[i]] + i
  }
  
  #classement des hauteurs
  for(i in 1:length(ordre_distRacineIN)){
    res[ordre_distRacineIN[i]] <- res[ordre_distRacineIN[i]] + i
  }
  
  classement <- sort(res, index.return= 1)$ix
  return (classement)
}




