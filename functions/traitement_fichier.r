root_dir <- getwd()
source(paste(root_dir, "/functions/protein_name.r", sep= ""))

list.of.packages <- c("seqinr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/")

library(seqinr)

##############################################
##  LECTURE ET CHARGEMENT MATRICE DISTANCES ##
##############################################

#' Remplie la partie superieure de la matrice avec sa partie inferieure, la diagonale est mise a 0
#'
#' @param matriceData : matrice de distances
#'
#' @return matrice de distance complete
#' @export Remplit la partie superieure de la matrice avec sa partie inferieure, la diagonale est mise a 0
#'
#' @examples data_403 <- remplacementNApar0(data_403)
remplacementNApar0 <- function(matriceData){
  tailleMat <- length(matriceData[1,])
  for(i in 1:tailleMat){
    for(j in 1:tailleMat){
      if(j > i) matriceData[i,j] <- matriceData[j,i]
      else if(i == j) matriceData[j,i] <- 0
    }
  }
  return (matriceData)
}

#' Place les donnees du dataframe dans une matrice de distances sans les noms des proteines
#'
#' @param fichier : nom du fichier .raw a lire
#' @param nb_sequence : nombre de sequences dans le fichier
#'
#' @return matrice de distances issue de la matrice diagonale inferieure du fichier (matrice de reels)
#' @export Place les donnees du dataframe dans une matrice de distances sans les noms des proteines
#'
#' @examples data_403 <- chargement_fichier("data/403_VLD_dist.raw", 403)
chargement_fichier <- function(fichier, nb_sequence){
  mydata <- read.table(fichier, sep="\t", dec= '.', fill= TRUE)
  
  data <- matrix(nrow= nb_sequence, ncol= nb_sequence)
  for(i in 1:nb_sequence){
    #on ignore toujours la premiere colonne
    for(j in 2:nb_sequence){
      if(2 == i){
        #ligne 2 = premiere vraie ligne
        data[i-1,j-1] <- 0
      } else if(1 == i){
        #premiere ligne = derniere vraie ligne
        data[nb_sequence,j-1] <- mydata[i,j]
      } else {
        #premiere ligne a val > 0 - 3e ligne
        data[i-1,j-1] <- mydata[i,j]
      }
    }
  }
  data <- remplacementNApar0(data)
  return (data)
}

#' Recupere les noms des proteines dans le fichier passe en parametre
#'
#' @param fichier : nom du fichier a lire (fichier .raw)
#'
#' @return noms des proteines du fichier (vecteur de chaines de caracteres)
#' @export Recupere les noms des proteines dans le fichier passe en parametre
#'
#' @examples names403 <- getProteinNames("data/403_VLD_dist.raw")
getProteinNames <- function(fichier){
  mydata <- read.table(fichier, sep="\t", dec= '.', fill= TRUE)
  names <- c()
  #on ajoutera la "premiere proteine" a la fin car c'est la derniere en fait
  for(i in 2:length(mydata)-1){
    names <- append(names, as.character(mydata[i,1]))
  }
  names <- append(names, as.character(mydata[1,1]))
  return (names)
}

#' Ecrit dans des fichiers propres aux groupes les noms des proteines qui les composent
#'
#' @param groupes : liste des groupes (pam$clustering)
#' @param noms_prot : noms des proteines
#' @param rd : ordre dans lequel les proteines ont ete selectionnees dans la liste de 1751
#'
#' @return RIEN - cree les fichiers .txt correspondants aux groupes retournes par pam
#' @export Ecrit dans des fichiers propres aux groupes les noms des proteines qui les composent
#'
#' @examples ecrire_fic_groupe(myPam$clustering, names1751, rd)
ecrire_fic_groupe <- function(groupes, noms_prot, rd){
  nb_different_groupes <- getNbGroups(groupes)
  for(i in 1:nb_different_groupes){
    nom_fic <- paste("protGr", i, sep= "")
    rep <- paste("resultats/", nom_fic, ".txt", sep= "")
    sink(rep, append=TRUE)
    for(j in 1:length(groupes)){
      if(i == groupes[j]){
        print(noms_prot[rd[j]])
      }
    }
    sink(NULL)
  }
}

#' Ecrit dans des fichiers les noms des proteines triees par groupes par le pam
#'
#' @param groupes : resultat de pam$clustering (vecteur d'entiers)
#' @param noms_prot : nombre de proteines (entier)
#'
#' @return RIEN - cree les fichiers .txt correspondants aux groupes retournes par pam
#' @export Ecrit dans des fichiers les noms des proteines triees par groupes par le pam
#'
#' @examples ecrire_fic_clustering(myPam$clustering, 403)
ecrire_fic_clustering <- function(groupes, noms_prot){
  nb_different_groupes <- getNbGroups(groupes)
  print(length(groupes))
  print(nb_different_groupes)
  for(i in 1:nb_different_groupes){
    nom_fic <- paste("protGr", i, sep="")
    rep <- paste("resultats/clustering/",nom_fic, ".txt", sep="")
    sink(rep, append=FALSE)
    for(j in 1:length(groupes)){
      if(i == groupes[j]){
        print(noms_prot[j])
      }
    }
    sink(NULL)
  }
}

#' Sauvegarde l'image d'une heatmap et la stocke dans le dossier '/resultats/heatmaps'
#'
#' @param heatmap : plot de la heatmap a enregistrer en .png
#' @param taille_hM : identifiant de la heatmap (entier)
#'
#' @return RIEN - cree un fichier .png
#' @export Sauvegarde l'image d'une heatmap et la stocke dans le dossier '/resultats/heatmaps'
#'
#' @examples ajouter_fic_heatmap(hM, 403)
ajouter_fic_heatmap <- function(heatmap, taille_hM){
  currentDir <- getwd()
  setwd("~/R/resultats/heatmaps/")
  
  x <- "ls | wc -l"
  nb_hM <- as.numeric(system(x, intern = TRUE)) + 1 # nombre de heatmaps deja stockees + 1
  new_heatmap <- paste("hM", taille_hM, "_", nb_hM, sep="")
  new_image <- paste(new_heatmap, ".png", sep="")
  
  if(!file.exists(new_image)){
    png(filename=new_image, width=500, height=500)
    heatmap
    dev.off()
  }
  
  setwd(currentDir)
}

#' Cree et stocke une heatmap (fichier .png) de 'matrice' dans le dossier 'resultats/heatmaps'
#'
#' @param dir_name : dossier racine du programme
#' @param file_name : nom du fichier a creer
#' @param matrice : cible de la heatmap
#'
#' @return heatmap (fichier .png) de 'matrice' dans le dossier 'resultats/heatmaps'
#' @export Cree et stocke une heatmap (fichier .png) de 'matrice' dans le dossier 'resultats/heatmaps'
#'
#' @examples ecriture_fichier_heatmap("~/home/user/", "monfichierHeatmap.png", matRob)
ecriture_fichier_heatmap <- function(dir_name, file_name, matrice){
  dirName <- paste(dir_name, "/resultats/heatmaps/", sep= "")
  if(!dir.exists(dirName)) dir.create(dirName)
  
  fileName <- paste(dirName, file_name, sep= "")
  if(!file.exists(fileName)){
    png(filename=fileName, width=500, height=500)
    heatmap(matrice)
    dev.off()
  }
}

##################################################
############  FONCTIONS CLIQUE ###################
##################################################


#' Creation du fichier de sommets utilise pour la clique
#'
#' @param root_dir: emplacement racine du programme
#' @param fic_sommets : nom du fichier 'sommets' souhaite (chaine de caracteres)
#' @param names : noms des proteines a stocker (vecteur de chaines de caracteres)
#'
#' @return RIEN - cree un fichier .txt
#' @export Creation du fichier de sommets utilise pour la clique
#'
#' @examples creer_fichier_sommets("sommets_file.txt", names403)
creer_fichier_sommets <- function(root_dir, fic_sommets, names){
  currentDir <- getwd()
  record_dir <- paste(root_dir, "/data/clique/", sep= "")
  setwd(record_dir)
  
  sink(fic_sommets, append= FALSE)
  cat("nom_prot regne type milieu\n")
  for(name in names){
    regne <- get_regne(name)
    type <- get_type(name)
    milieu <- get_milieu(name)
    cat(paste(name, regne, type, milieu), "\n")
  }
  sink(NULL)
  
  setwd(currentDir)
}

#' Creation du fichier d'arcs utilise pour la creation de clique
#'
#' @param root_dir: emplacement racine du programme
#' @param fic_arcs : nom du fichier 'arcs' a creer (chaine de caracteres)
#' @param matRobustesse : matrice de robustesse (matrice d'entiers)
#' @param names : noms des proteines (vecteur de chaines de caracteres)
#' @param force : lien qui doit lier les individus entre eux
#'
#' @return RIEN - cree un fichier .txt
#' @export Creation du fichier d'arcs utilise pour la creation de clique
#'
#' @examples creer_fichier_arcs("arc_file.txt", matRob, names403)
creer_fichier_arcs <- function(root_dir, fic_arcs, matRobustesse, names, force){
  currentDir <- getwd()
  data_dir <- paste(root_dir, "/data/clique/", sep= "")
  setwd(data_dir)
  
  sink(fic_arcs, append= FALSE)
  cat("prot_1 prot_2 force\n")
  for(i in 1:length(matRobustesse[1,])){
    for(j in 1:i){
      if((i != j) && (force <= matRobustesse[i,j])){
        prot_i <- names[i]
        prot_j <- names[j]
        robustesse <- matRobustesse[i,j]
        cat(paste(prot_i, prot_j, robustesse, sep=" "), "\n")
      }
    }
  }
  sink(NULL)
  setwd(currentDir)
}


#' Sauvegarde l'image d'une clique et la stocke dans le dossier clique
#'
#' @param clique_name : nom de fichier desire
#' @param clique : resultat de la fonction plot
#'
#' @return RIEN - cree un fichier .png ()
#' @export Sauvegarde l'image d'une clique et la stocke dans le dossier clique
#'
#' @examples ajouter_fic_clique("clique403.png", plot(clique))
ajouter_fic_clique <- function(clique_name, clique){
  currentDir <- getwd()
  setwd("~/R/resultats/clique/")
  
  new_image <- clique_name
  if(file.exists(new_image)){
    new_image <- paste(new_image, "(1)", sep="")
    count <- 2
    while(file.exists(new_image)){
      new_image <- paste(clique_name, "(", count, ")", sep="")
      count <- count+1
    }
  }
  png(filename=new_image, width=500, height=500)
  clique
  dev.off()
  
  setwd(currentDir)
}

##################################################
##  CORRECTION DU FICHIER DES NOMS DE PROTEINES ##
##################################################

#' Retourne une matrice (nb_prot x 2) avec sur chaque ligne l'ancien nom et l'eventuel nouveau nom
#'
#' @param fic_name: nom du fichier de correction (chaine de caracteres)
#' @param nb_prot : nombre de proteines dans le fichier (entier)
#'
#' @return matrice (nb_prot x 2) avec sur chaque ligne l'ancien nom et l'eventuel nouveau nom
#' @export Retourne une matrice (nb_prot x 2) avec sur chaque ligne l'ancien nom et l'eventuel nouveau nom
#'
#' @examples tab_real_names <- get_tab_real_names("new_names1751.txt", 1751)
get_tab_real_names <- function(fic_name, nb_prot){
  currentDir <- getwd()
  setwd("~/R/data/")
  
  frame_noms <- read.table(fic_name, header= FALSE, sep="\n", comment.char= "")
  noms <- matrix(nrow= nb_prot, ncol= 2)
  for(i in 1:nb_prot+1) #parcours des lignes du fichier lu par read.table
  {
    ligne <- unlist(strsplit(as.character(frame_noms[i,]), "\t")) #traite un vector et pas une liste avec 'unlist'
    noms[i-1, 1] <- ligne[1] #ancien nom de proteine
    if(2 == length(ligne)) noms[i-1, 2] <- ligne[2]
    else noms[i-1, 2] <- "" #pas de nouveau nom de proteine donc chaine vide
  }
  setwd(currentDir)
  return (noms)
}


#' Verifie que le fichier de correction des noms de proteines est correct selon 3 criteres :
#'  1) ancien nom different de nouveau nom
#'  2) base du nom de la proteine faisant bien 6 caracteres (AcamarTOP... est correct, nchar("Acamar") == 6 )
#'  3) aucune proteine dupliquee parmi la liste
#'
#' @param tab_real_names : resultat de 'get_tab_real_names'
#'
#' @return ensemble des noms des proteines (vecteur de chaines de caracteres)
#' @export Verifie que le fichier de correction des noms de proteines est correct selon certains criteres
#'
#' @examples real_names <- get_real_names(tab_real_names)
get_real_names <- function(tab_real_names){
  nb_egalites <- 0
  nb_mauvais_format <- 0
  nb_duplication <- 0
  
  real_names <- c()
  taille <- length(tab_real_names[,1])
  for(i in 1:taille){
    old <- tab_real_names[i, 2]
    new_diez <- tab_real_names[i, 1]
    new <- substr(new_diez, 2, nchar(new_diez))
    
    if("" == old) # Pas de changement de nom
    {
      new <- tab_real_names[i, 1]
      real_names <- append(real_names, new)
    } else # Changement de nom
    {
      if(old == new) {
        print(paste("Ligne ", i, " => ancien mot = nouveau mot"))
        nb_egalites <- nb_egalites + 1
      }
      if(!hasGoodNameFormat(new)){
        print(paste("Ligne ", i, " => ", new, " ne fait pas 6 caracteres"))
        nb_mauvais_format <- nb_mauvais_format + 1
      }
      if(new %in% real_names){
        print(paste("Ligne ", i, " => proteine deja presente ailleurs"))
        nb_duplication <- nb_duplication + 1
      }
      else real_names <- append(real_names, new)
    }
  }
  
  print("")
  print(paste("Nombre d'egalites : ", nb_egalites))
  print(paste("Nombre de mauvais formats : ", nb_mauvais_format))
  print(paste("Nombre duplications : ", nb_duplication))
  
  # Ligne 129 du fichier texte modifiee : #ArtE1TOPz_M_B ==> #ArtE1.TOPz_M_B
  
  return (real_names)
}


#' Retourne le prochain 'fic_name' non-existant dans le dossier '~R/resultats/clique/' se basant sur 'fic_name'
#'
#' @param fic_name : nom du fichier voulu
#'
#' @return fic_name si le fichier n'existe pas, fic_name(X) si le fichier existe (avec X un entier)
#' @export Retourne le prochain 'fic_name' non-existant dans le dossier '~R/resultats/clique/' se basant sur 'fic_name'
#'
#' @examples fic_name <- get_next_filename("clique.txt") ; si "clique.txt" existe, "clique(1).txt" est retourne
get_next_filename <- function(fic_name){
  currentDir <- getwd()
  setwd("~/R/resultats/clique/")
  
  new_image <- fic_name
  nom_fichier <- substr(new_image, 1, nchar(new_image)-4) #nom de fichier sans l'extension .png
  extension <- substr(new_image, nchar(new_image)-3, nchar(new_image)) #extension .png ou .jpg
  
  if(file.exists(new_image)){
    new_image <- paste(nom_fichier, "(1)", extension, sep="")
    count <- 2
    while(file.exists(new_image)){
      new_image <- paste(nom_fichier, "(", count, ")", extension, sep="")
      count <- count+1
    }
  }
  setwd(currentDir)
  return (new_image)
}

#############################################
#######  ECRITURE FICHIERS DE GROUPES #######
#############################################

#' Cree un fichier classant les proteines en differents groupes
#'
#' @param base_dir : repertoire racine du programme
#' @param folder_name : nom du dossier ou se placer (chaine de caracteres)
#' @param fic_name : nom du fichier desire (chaine de caracteres)
#' @param amisX : ensemble des groupes d'amis (liste de liste d'entiers)
#' @param namesX : noms des proteines presentes dans les groupes d'amis (vecteur de chaines de caracteres)
#'
#' @return RIEN - cree des fichiers .txt
#' @export Cree un fichier classant les proteines en differents groupes
#'
#' @examples ecriture_fichier_groupes("coupe", "coupe10.txt", amis403, names403)
ecriture_fichier_groupes <- function(base_dir, folder_name, fic_name, amisX, namesX){
  currentDir <- getwd()
  baseDir <- paste(base_dir, "/resultats/", sep= "")
  setwd(baseDir)
  
  #test de l'existence du dossier dans le repertoire 'baseDir'
  fullDirPath <- paste(baseDir, folder_name, sep= "")
  testDir <- paste("[ -d", fullDirPath,"]")
  existDir <- system(testDir)
  
  #creation du dossier s'il n'existe pas
  if(existDir >= 1){
    newDir <- paste("mkdir", folder_name)
    system(newDir)
  }
  
  #positionnement dans le dossier
  myDir <- paste(baseDir, folder_name, sep= "")
  setwd(myDir)
  
  nb_groupes <- length(amisX)
  nb_prot <- length(namesX)
  
  #creation du fichier
  sink(fic_name, append=FALSE)
  cat("Nombre de groupes : ", nb_groupes, "\n")
  cat("Ligne : Groupe ID_proteine nom_proteine\n")
  
  #groupe 0 - groupe particulier n'ayant pas d'amis
  singletons <- find_singletons(amisX, nb_prot)
  cat("##############################################\n")
  cat("# Singletons - ", length(singletons), " proteines\n")
  for(s in 1:length(singletons)){
    cat("0 ", singletons[s], namesX[singletons[s]],"\n")
  }
  
  #groupes d'amis
  for(i in 1:nb_groupes){
    local_friends <- sort(unlist(amisX[[i]]))
    taille_groupe <- length(local_friends)
    cat("##############################################\n")
    cat("# Groupe ", i, " - ", taille_groupe, " proteines\n")
    for(j in 1:taille_groupe){
      local_friend <- local_friends[j]
      cat(i, " ", local_friend, namesX[local_friend],"\n")
    }
  }
  sink(NULL)
  setwd(currentDir)
}

#' Applique hclust a la matrice, puis cutree et enregistre les groupes trouves dans le dossier voulu
#'
#' @param mat : matrice de distance ou de robustesse a utiliser (matrice de reels ou d'entiers)
#' @param typeMat : type de matrice (distance, robustesse etc) (utilise pour le nom du fichier cree) (chaine de caracteres)
#' @param methode : nom de la methode a utiliser (ward.D, ward.D2, single, complete, average, mcquitty, median, centroid) "ward.D2" par defaut (chaine de caracteres)
#' @param coupe : nombre de groupes a former avec hclust (entier)
#' @param folder_name : nom du dossier dans lequel enregistrer les groupes formes (chaine de caracteres)
#' @param namesX : noms des proteines (vecteur de chaines de caracteres)
#'
#' @return resultat du cutree (vecteur d'entiers)
#' @export Applique hclust a la matrice, puis cutree et enregistre les groupes trouves dans le dossier voulu
#'
#' @examples groupes <- cutAndWrite(matRob403, "robustesse", "ward.D", 10, "~/home/user/", "coupe", names403)
cutAndWrite <- function(mat, typeMat, methode= "ward.D2", nb_groupes, base_dir, folder_name, namesX){
  retMat <- mat
  clust_retMat <- hclust(as.dist(retMat), method = methode, members= NULL)
  plot_name <- paste(typeMat, length(mat[1,]), "_", methode, sep= "")
  plot(clust_retMat, main = plot_name) #affichage du dendogramme associe
  
  print("")
  print("L'arbre a couper se trouve dans le dossier resultats")
  print("Combien de groupes voulez-vous ?")
  nb_groupes <- read(nmax=1)
  
  coupe <- cutree(clust_retMat, k= nb_groupes)
  amis_coupe <- build_friend_list(coupe, nb_groupes)
  
  #ecriture des fichiers de groupes formes par cutree
  nb_prot <- length(namesX)
  nom_fichier <- paste(typeMat, nb_prot, "_", methode, "_coupe", nb_groupes, ".txt", sep= "")
  ecriture_fichier_groupes(base_dir, folder_name, nom_fichier, amis_coupe, namesX)
  
  return (coupe)
}

##################################################
##########  FONCTIONS OCCURRENCE #################
##################################################

#' Extrait les differents caracteres du fichier .fasta lu par 'read.fasta'
#'
#' @param struct_fasta : resultat de la fonction read.fasta
#'
#' @return alphabet du fichier fasta lu (vecteur de chaines de caracteres)
#' @export Extrait les differents caracteres du fichier .fasta lu par 'read.fasta'
#'
#' @examples alphabet <- getRealAlphabet(struct_fasta)
getRealAlphabet <- function(struct_fasta){
  alphabet <- list()
  
  for(i in 1:length(struct_fasta)){
    identifiants <- strsplit(struct_fasta[[i]], split= " ")
    
    for(j in 1:length(identifiants)){
      alphabet <- union(alphabet, identifiants[[j]])
    }
  }
  alphabet_vector <- unlist(alphabet)
  return (alphabet_vector)
}

#' Retourne la position du caractere dans l'alphabet
#'
#' @param alphabet : ensemble des caracteres (liste de chaines de caracteres)
#' @param caract : caractere dont on veut l'indice (chaine de caractere)
#'
#' @return indice du caractere dans l'alphabet (entier)
#' @export Retourne la position du caractere dans l'alphabet
#'
#' @examples positionE10 <- getCaractLocation(alphabet, "E10")
getCaractLocation <- function(alphabet, caract){
  for(i in 1:length(alphabet)){
    if(alphabet[i] == caract) return (i)
  }
  return (0)
}

#' Retourne les noms et le nombre d'occurrence de chaque caractere du fichier .fasta
#'
#' @param file_name : nom du fichier .fasta (doit se trouver dans '~/R/data/') (chaine de caracteres)
#'
#' @return caracteres presents dans 'file_name' + occurrences (liste : vecteur de chaines de caracteres + vecteur d'entiers)
#' @export Retourne le nombre d'occurrence de chaque caractere du fichier .fasta
#'
#' @examples alphaOcc <- getOccCaract("myfile.fasta") ; alphabet <- alphaOcc[[1]] ; occurrence <- alphaOcc[[2]]
getOccCaract <- function(file_path){
#   baseDir <- "~/R/data/"
  # file_path <- paste(baseDir, file_name, sep= "") #chemin complet du fichier a lire
  struct_fasta <- read.fasta(file_path, as.string = TRUE, seqonly = TRUE, forceDNAtolower= FALSE) #ensemble des sequences fasta
  taille_fasta <- length(struct_fasta) #nombre de sequences fasta
  
  alphabet <- getRealAlphabet(struct_fasta)
  taille_alphabet <- length(alphabet)
  occurrence <- c(rep(0, taille_alphabet)) #vecteur qui compte le nombre d'occurrence de chaque elem de l'alphabet
  
  for(j in 1:taille_fasta) #pour toutes les sequences de struct_fasta
  {
    print(paste("recuperation occurrence", j, "/", taille_fasta))
    identifiants <- unlist(strsplit(struct_fasta[[j]], split= " "))
    for(k in 1:length(identifiants)) #pour tous les elements de la sequence j
    {
      if(0 == getCaractLocation(alphabet, identifiants[k])) print("caractere non trouve")
      occurrence[getCaractLocation(alphabet, identifiants[k])] <- occurrence[getCaractLocation(alphabet, identifiants[k])] + 1
    }
  }
  
  occCaract <- list(alphabet, occurrence)
  return (occCaract)
}

#' Ecrit dans un fichier les differents caracteres d'un alphabet et leurs occurrences respectives
#'
#' @param fic_name : nom du fichier desire (sera cree dans '~/R/resultats/occurrence/') (chaine de caracteres)
#' @param occCaract : resultat de la fonction 'getOccCaract' (liste : vecteur de chaines de caracteres + vecteur d'entiers)
#'
#' @return RIEN - cree un fichier .txt
#' @export Ecrit dans un fichier les differents caracteres d'un alphabet et leurs occurrences respectives
#'
#' @examples ecriture_fichier_occurrence("fichierOcc.txt", alphaOcc)
ecriture_fichier_occurrence <- function(root_dir, fic_name, occCaract){
  currentDir <- getwd()
  path <- paste(root_dir, "/resultats/occurrence/", sep= "")
  setwd(path)
  
  sink(fic_name, append= FALSE)
  
  for(boucle in 1:length(occCaract[[1]])){
    cat(occCaract[[1]][boucle],"-", occCaract[[2]][boucle], "\n", sep= "")
  }
  
  sink(NULL)
  setwd(currentDir)
}

#' Recupere l'alphabet et l'occurrence de chaque caractere, les trie et les ecrit dans un fichier
#'
#' @param fasta_file : fichier fasta de donnees
#' @param new_file_name : nom du fichier 'alphabet et occurrences' a creer
#'
#' @return l'alphabet et l'occurrence de chaque caractere
#' @export Recupere l'alphabet et l'occurrence de chaque caractere, les trie et les ecrit dans un fichier
#'
#' @examples occCaractTest <- getSortAndWriteOccCaract("403_seq.vld.fasta", "occCaract403")
getSortAndWriteOccCaract <- function(fasta_file, new_file_name){
  occCaract <- getOccCaract(fasta_file) #recuperation de l'alphabet + occurrences
  new_occCaract <- triAlphabetDecroissant(occCaract) #tri alphabet et occurrences
  ecriture_fichier_occurrence(new_file_name, new_occCaract) #ecriture fichier occurrences
  
  return (occCaract)
}

getSubAndRetMat <- function(mat, nb_prot, percent){
  remaining_ones <- sort(sample(1:nb_prot, nb_prot - (nb_prot*percent/100)))
  subMat <- retournementMat(mat[remaining_ones, remaining_ones], mat[1,1])
  
  return (subMat)
}

#' Tri de l'alphabet sur la taille des caracteres
#'
#' @param occCaract : alphabet et occurrence
#'
#' @return alphabet et occurrence tries sur la taille des caracteres
#' @export Tri de l'alphabet sur la taille des caracteres
#'
#' @examples occCaract <- triAlphabet(occCaract)
triAlphabet <- function(occCaract){
  alphabet <- occCaract[[1]]
  occurrence <- occCaract[[2]]
  
  new_alphabet <- c()
  new_occurrence <- c()
  
  for(i in 1:length(alphabet)){
    if(1 == nchar(alphabet[i])) #tous les caracteres non transformes sont ajoutes au debut
    {
      new_alphabet <- append(new_alphabet, alphabet[i])
      new_occurrence <- append(new_occurrence, occurrence[i])
    }
  }
  
  for(i in 1:length(alphabet)){
    if(1 < nchar(alphabet[i])) #tous les caracteres transformes
    {
      new_alphabet <- append(new_alphabet, alphabet[i])
      new_occurrence <- append(new_occurrence, occurrence[i])
    }
  }
  
  res <- list(new_alphabet, new_occurrence)
  return (res)
}


#' Trie l'alphabet selon l'occurrence decroissante des caracteres
#'
#' @param occCaract : alphabet et occurrence (liste : vecteur de chaines de caracteres + vecteur d'entiers)
#'
#' @return alphabet trie par occurrence decroissante
#' @export Trie l'alphabet selon l'occurrence decroissante des caracteres
#'
#' @examples occCaract <- triAlphabetDecroissant(occCaract)
triAlphabetDecroissant <- function(occCaract){
  occDecroissante <- sort(unlist(occCaract[[2]]), decreasing= TRUE, index.return= 1)$ix
  
  alphabet <- unlist(occCaract[[1]])[occDecroissante]
  occurrence <- unlist(occCaract[[2]])[occDecroissante]
  
  res <- list(alphabet, occurrence)
  return (res)
}


#' Enregistre l'histogramme des occurrences dans un fichier .png
#'
#' @param root_dir : repertoire racine du programme
#' @param file_name : chemin du fichier a creer
#' @param occCaract : alphabet et occurrence (liste : vecteur de chaines de caracteres + vecteur d'entiers)
#'
#' @return RIEN - cree un fichier .png
#' @export Enregistre l'histogramme des occurrences dans un fichier .png
#'
#' @examples ecriture_fichier_hist("~/R/resultats/hist_occurrence403.png", occCaract)
ecriture_fichier_hist_occurrence <- function(root_dir, file_name, occCaract){
  currentDir <- getwd()
  path <- paste(root_dir, "/resultats/occurrence/", sep= "")
  setwd(path)
  
  occ <- occCaract[[2]] #occurrence des caracteres de l'alphabet
  filePath <- paste(path, "/", file_name, ".png", sep= "")
  png(filename= filePath)
  hist(table(data.frame(occ)), xlab= "Occurrence", main= paste("Histogramme de ", file_name, sep= ""))
  dev.off()
  
  setwd(currentDir)
}

##################################################
############  FONCTIONS FASTA ####################
##################################################

#' Recupere les noms et sequences des proteines presentes dans le fichier
#'
#' @param baseDir : repertoire racine du programme
#' @param file_name : chemin du fichier .fasta (doit se trouver dans "data/") (chaine de caracteres)
#'
#' @return noms et sequences des proteines (liste de deux vecteurs de chaines de caracteres)
#' @export Recuperation des noms et sequences proteiques d'un fichier fasta
#'
#' @examples nomsEtsequences <- getNamesAndSeq("~/home/user/", "myfile.fasta")
getNamesAndSeq <- function(baseDir, file_path){
  struct_fasta <- read.fasta(file_path, as.string = TRUE, forceDNAtolower= FALSE) #ensemble des sequences fasta
  
  taille <- length(struct_fasta)
  fullNames <- getAnnot(struct_fasta) #noms complets des proteines
  sequences <- vector(length= taille) #sequences completes
  for(i in 1:taille){
    fullNames[[i]] <- gsub("\t", " ", fullNames[[i]])
    sequences[i] <- struct_fasta[[i]][[1]]
  }
  res <- list(fullNames, sequences)
  return (res)
}


#' Retourne le vecteur des chaines creees en ayant decoupe 'chaine' en 'val' sous-chaines
#'
#' @param chaine : sequence a decouper (chaine de caracteres)
#' @param val : nombre de sous-chaines souhaite (entier)
#'
#' @return vecteur des chaines creees en ayant decoupe 'chaine' en 'val' sous-chaines
#' @export Retourne le vecteur des chaines creees en ayant decoupe 'chaine' en 'val' sous-chaines
#'
#' @examples ss-str <- wordCut("blablobla", 3) ; ss-str = ("bla", "blo", "bla")
wordCut <- function(chaine, val){
  taille <- (nchar(chaine)/val+1)
  res <- vector(length= taille)
  for(i in 0:taille-1){
    sous_chaine <- substr(chaine, i*val+1, i*val+val)
    res[i+1] <- sous_chaine
  }
  
  return (res)
}

#' Cree un fichier fasta avec les proteines de la liste 'amis'
#'
#' @param folder_name : nom du dossier cible
#' @param fic_name : nom du fichier fasta a creer
#' @param fullNamesAndSeq : resultat de getNamesAndSeq (liste de 2 vecteurs de chaines de caracteres)
#' @param amis : liste des proteines du groupe
#'
#' @return RIEN - cree un fichier fasta avec les proteines de la liste 'amis'
#' @export Cree un fichier fasta avec les proteines de la liste 'amis'
#'
#' @examples create_fichier_fasta("Fic403.fasta", fullNamesAndSeq_403, amis403)
create_fichier_fasta <- function(folder_name, fic_name, fullNamesAndSeq, amis){
  nb_prot <- length(amis)
  X <- 60 # valeur de decoupe des sequences pour l'affichage dans le fichier
  fic_path <- paste(folder_name, "/", fic_name, sep= "")
  
  sink(fic_path, append= FALSE)
  for(i in 1:nb_prot)
  {
    nom_prot <- fullNamesAndSeq[[1]][[amis[[i]]]]
    seq_prot <- fullNamesAndSeq[[2]][[amis[[i]]]]
    cat(nom_prot, "\n", sep= "")
    cut_seq <- wordCut(seq_prot, X) #decoupe la sequence en sous chaines de X caracteres
    for(j in 1:length(cut_seq))
    {
      sous_chaine <- cut_seq[j]
      if(sous_chaine != "") cat(sous_chaine, "\n", sep= "")
    }
  }
  sink(NULL)
}

#' Cree les fichiers fasta correspondants aux groupes d'amis dans differents dossiers
#'
#' @param folder_name : nom du dossier où enregistrer les fichiers fasta
#' @param fullNamesAndSeq : resultat de getNamesAndSeq (liste de 2 vecteurs de chaines de caracteres)
#' @param amis : ensemble des groupes d'amis trouves (liste de listes d'entiers)
#'
#' @return RIEN - cree des fichiers .fasta
#' @export Creation de fichiers .fasta lies aux groupes de proteines classees ensemble
#'
#'@examples ecriture_all_fichiers_fasta("fasta", nomsEtsequences, amis)
ecriture_all_fichiers_fasta <- function(folder_name, fullNamesAndSeq, amis){
  currentDir <- getwd()
  dirName <- folder_name
  if(!dir.exists(dirName)) dir.create(dirName)
  setwd(dirName)
  
  for(i in 1:length(amis))
  {
    fic_name <- paste("Groupe", i, ".fasta", sep= "")
    create_fichier_fasta(dirName, fic_name, fullNamesAndSeq, unlist(amis[[i]]))
  }
  
  setwd(currentDir)
}

##################################################
############  FONCTIONS VENN #####################
##################################################

#' Enregistre les png des diagrammes de Venn crees avec les 3 versions de coupe passees en parametre
#'
#' @param version_1 : premiere version de repartition des groupes
#' @param version_2 : deuxieme version de repartition des groupes
#' @param version_3 : troisieme version de repartition des groupes
#' @param nb_groupes : nombre de groupes crees pour la coupe
#'
#' @return png des diagrammes de Venn crees avec les 3 versions de coupe passees en parametre
#' @export Enregistre les png des diagrammes de Venn crees avec les 3 versions de coupe passees en parametre
#'
#' @examples ecriture_fichiers_venn(DATAcoupe10wardD2, RETcoupe10wardD2, SUBcoupe10wardD2, 10)
ecriture_fichiers_venn <- function(version_1, version_2, version_3, nb_groupes){
  
  noms_coupes <- c("DATA", "RET", "SUB")
  couleurs_coupes <- c("skyblue", "pink1", "mediumorchid")
  
  for(i in 1:nb_groupes) #differents groupes de la coupe
  {
    saveDir <- "~/R/resultats/venn/"
    if(!dir.exists(saveDir)) dir.create(saveDir)
    
    file_name <- paste("Venn_groupe_", i, ".png", sep= "")
    file_path <- paste(saveDir, file_name, sep= "")
    
    D1 <- getGroupeCoupe(version_1, i)
    R1 <- getGroupeCoupe(version_2, i)
    S1 <- getGroupeCoupe(version_3, i)
    
    taille_D1 <- length(D1)
    taille_R1 <- length(R1)
    taille_S1 <- length(S1)
    
    intersectionDR <- length(intersect(D1, R1))
    intersectionDS <- length(intersect(D1, S1))
    intersectionSR <- length(intersect(S1, R1))
    
    intersectionDRS <- length(intersect(intersect(D1, R1), S1))
    
    #enregistrement du diagramme de Venn
    png(file_path)
    grid.newpage()
    draw.triple.venn(area1= taille_D1, area2= taille_R1, area3= taille_S1,
                     n12= intersectionDR, n13= intersectionDS, n23= intersectionSR,
                     n123= intersectionDRS,
                     category= noms_coupes,
                     fill= couleurs_coupes)
    
    dev.off()
  }
}