##############################################
##  LECTURE ET CHARGEMENT MATRICE DISTANCES ##
##############################################

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
  return (data)
}

#' Remplit la partie superieure de la matrice avec sa partie inferieure, la diagonale est mise a 0
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
  for(i in 2:length(mydata)){
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
    nom_fic <- paste("protGr", i, sep="")
    rep <- paste("resultats/",nom_fic, ".txt", sep="")
    sink(rep, append=TRUE)
    for(j in 1:length(groupes)){
      if(i == groupes[j]){
        print(names1751[rd[j]])
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

##################################################
############  FONCTIONS CLIQUE ###################
##################################################


#' Creation du fichier de sommets utilise pour la clique
#'
#' @param fic_sommets : nom du fichier 'sommets' souhaite (chaine de caracteres)
#' @param names : noms des proteines a stocker (vecteur de chaines de caracteres)
#'
#' @return RIEN - cree un fichier .txt
#' @export Creation du fichier de sommets utilise pour la clique
#'
#' @examples creer_fichier_sommets("sommets_file.txt", names403)
creer_fichier_sommets <- function(fic_sommets, names){
  currentDir <- getwd()
  setwd("~/R/data/clique/")
  
  sink(fic_sommets, append= FALSE)
  cat("nom_prot type part milieu\n")
  for(name in names){
    type <- get_type(name)
    part <- get_particularite(name)
    milieu <- get_milieu(name)
    cat(sprintf(paste(name, type, part, milieu), sep=" "), "\n")
  }
  sink(NULL)
  
  setwd(currentDir)
}

#' Creation du fichier d'arcs utilise pour la creation de clique
#'
#' @param fic_arcs : nom du fichier 'arcs' a creer (chaine de caracteres)
#' @param matRobustesse : matrice de robustesse (matrice d'entiers)
#' @param names : noms des proteines (vecteur de chaines de caracteres)
#'
#' @return RIEN - cree un fichier .txt
#' @export Creation du fichier d'arcs utilise pour la creation de clique
#'
#' @examples creer_fichier_arcs("arc_file.txt", matRob, names403)
creer_fichier_arcs <- function(fic_arcs, matRobustesse, names){
  currentDir <- getwd()
  setwd("~/R/data/clique/")
  
  sink(fic_arcs, append= FALSE)
  cat("prot_1 prot_2 force\n")
  for(i in 1:length(matRobustesse[1,])){
    for(j in 1:i){
      if(i != j){
        prot_i <- names[i]
        prot_j <- names[j]
        robustesse <- matRobustesse[i,j]
        cat(sprintf(paste(prot_i, prot_j, robustesse, sep=" ")), "\n")
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
#' @param folder_name : nom du dossier ou se placer (chaine de caracteres)
#' @param fic_name : nom du fichier desire (chaine de caracteres)
#' @param amisX : ensemble des groupes d'amis (liste de liste d'entiers)
#' @param namesX : noms des proteines presentes dans les groupes d'amis (vecteur de chaines de caracteres)
#'
#' @return RIEN - cree des fichiers .txt
#' @export Cree un fichier classant les proteines en differents groupes
#'
#' @examples ecriture_fichier_groupes("coupe", "coupe10.txt", amis403, names403)
ecriture_fichier_groupes <- function(folder_name, fic_name, amisX, namesX){
  currentDir <- getwd()
  baseDir <- "~/R/resultats/"
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
#' @return caracteres presents dans 'file_name' + occurrences (liste : liste de chaines de caracteres + liste d'entiers)
#' @export Retourne le nombre d'occurrence de chaque caractere du fichier .fasta
#'
#' @examples alphaOcc <- getOccCaract("myfile.fasta") ; alphabet <- alphaOcc[[1]] ; occurrence <- alphaOcc[[2]]
getOccCaract <- function(file_name){
  baseDir <- "~/R/data/"
  file_path <- paste(baseDir, file_name, sep= "") #chemin complet du fichier a lire
  struct_fasta <- read.fasta(file_path, as.string = TRUE, seqonly = TRUE, forceDNAtolower= FALSE) #ensemble des sequences fasta
  taille_fasta <- length(struct_fasta) #nombre de sequences fasta
  
  alphabet <- getRealAlphabet(struct_fasta)
  taille_alphabet <- length(alphabet)
  occurrence <- c(rep(0, taille_alphabet)) #vecteur qui compte le nombre d'occurrence de chaque elem de l'alphabet
  
  for(j in 1:taille_fasta) #pour toutes les sequences de struct_fasta
  {
    print(paste("sequence", j, "/", taille_fasta))
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
#' @param occCaract : resultat de la fonction 'getOccCaract' (liste : liste de chaines de caracteres + liste d'entiers)
#'
#' @return RIEN - cree un fichier .txt
#' @export Ecrit dans un fichier les differents caracteres d'un alphabet et leurs occurrences respectives
#'
#' @examples ecriture_fichier_occurrence("fichierOcc.txt", alphaOcc)
ecriture_fichier_occurrence <- function(fic_name, occCaract){
  currentDir <- getwd()
  setwd("~/R/resultats/occurrence/")
  
  sink(fic_name, append= FALSE)
  
  for(boucle in 1:length(occCaract[[1]])){
    cat(occCaract[[1]][boucle],"-", occCaract[[2]][boucle], "\n", sep= "")
  }
  
  sink(NULL)
  setwd(currentDir)
}

#' Recupere les noms et sequences des proteines presentes dans le fichier
#'
#' @param file_name : nom du fichier .fasta (doit se trouver dans ~/R/data/) (chaine de caracteres)
#'
#' @return noms et sequences des proteines (liste de deux vecteurs de chaines de caracteres)
#' @export Recuperation des noms et sequences proteiques d'un fichier fasta
#'
#' @examples nomsEtsequences <- getNamesAndSeq("myfile.fasta")
getNamesAndSeq <- function(file_name){
  baseDir <- "~/R/data/"
  file_path <- paste(baseDir, file_name, sep= "") #chemin complet du fichier a lire
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


#' Cree les fichiers fasta correspondants aux groupes d'amis dans differents dossiers
#'
#' @param fullNamesAndSeq : resultat de getNamesAndSeq (liste de 2 vecteurs de chaines de caracteres)
#' @param amis : ensemble des groupes d'amis trouves (liste de listes d'entiers)
#'
#' @return RIEN - cree des fichiers .fasta
#' @export Creation de fichiers .fasta lies aux groupes de proteines classees ensemble
#'
#'@examples ecriture_all_fichiers_fasta(nomsEtsequences, amis)
ecriture_all_fichiers_fasta <- function(fullNamesAndSeq, amis){
  
  dirName <- "~/R/resultats/coupe/coupe10"
  if(!dir.exists(dirName)) dir.create(dirName)
  
  
  
  
}

