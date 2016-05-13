#lit le fichier .raw avec la derniere proteine remontee en premiere position
#et avec le nom de la proteine en premiere position
lecture_raw <- function(mydata, nb_sequence){
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

#place les donnees du dataframe dans une matrice de distances sans les noms des proteines
chargement_fichier <- function(fichier, nb_seq){
  mydata <- read.table(fichier, sep="\t", dec= '.', fill= TRUE)
  data <- lecture_raw(mydata, nb_seq)
  return (data)
}

#remplace les valeurs NA par des 0 dnas la matrice pour faciliter les calculs
remplacementNApar0 <- function(matriceData){
  tailleMat <- length(matriceData[1,])
  for(i in 1:tailleMat){
    for(j in 1:tailleMat){
      if(j > i) matriceData[i,j] <- matriceData[j,i]
      # else if(is.na(matriceData[i,j])) matriceData[i,j] <- 0
      else if(i == j) matriceData[j,i] <- 0
    }
  }
  return (matriceData)
}

#recupere les noms des proteines dans le fichier qu'on passe en parametre
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

#ecrit dans des fichiers propres aux groupes les noms des proteines qui les composent
#groupes = liste des groupes (pam$clustering); 
#noms_prot = noms des proteines ;
#rd = ordre dans lequel les proteines ont ete selectionnees dans la liste de 1751
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

#ecrit dans des fichiers les noms des proteines triees par groupes par le pam
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

#sauvegarde l'image d'une heatmap et la stocke dans le dossier heatmaps
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


# Creation du fichier de sommets utilise pour la clique
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

# Creation du fichier d'arcs utilise pour la clique
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


#sauvegarde l'image d'une heatmap et la stocke dans le dossier heatmaps
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


get_next_filename <- function(fic_name){
  currentDir <- getwd()
  setwd("~/R/resultats/clique/")
  
  new_image <- fic_name
  nom_fichier <- substr(new_image, 1, nchar(new_image)-4)
  extension <- substr(new_image, nchar(new_image)-3, nchar(new_image))
  
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

build_fasta_tab <- function(file_name){
  baseDir <- "~/R/data/"
  file_path <- paste(baseDir, file_name, sep= "")
  
  struct_fasta <- read.fasta(file_path)
  #   print(names(fasta_fic)[1])
  #   print(fasta_fic[[1]][[1]])
  taille_fasta <- length(struct_fasta)
  
  for(elem in struct_fasta){
    print(elem)
  }
  
  alphabet <- getAlphabet(struct_fasta)
  taille_alphabet <- length(alphabet)
  occurrence <- c(rep(0, taille_alphabet))
  
  print(length(occurrence))
  
  print(unlist(alphabet))
  print(paste("taille alphabet = ", taille_alphabet))
  
  for(i in 1:taille_fasta) #pour chaque sequence fasta
  {
    for(j in 1:length(struct_fasta[[i]]))#pour chaque caractere de cette sequence
    {
      caract <- struct_fasta[[i]][[j]]
      if(" " != caract){
        print(paste("taille sequence ", length(struct_fasta[[i]])))
        char_checked= FALSE
        k <- 1
        while(!char_checked){
          print(paste("k=", k))
          print(paste("struct=", struct_fasta[[i]][[j]]))
          print(paste("alpha=", alphabet[[k]]))
          
          
          if(alphabet[[k]] == caract) #caractere trouve
          {
            occurrence[k] <- occurrence[k] + 1
            char_checked= TRUE
          } else {
            k <- k + 1
          }
        }
      }
    }
  }
  
  #calcul des occurrences des caracteres de l'alphabet
  
  res_mat <- list(alphabet, as.list(occurrence))
  
  return (res_mat)
}

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
  
  sink(fic_name, append=FALSE)
  cat("Nombre de groupes : ", nb_groupes, "\n")
  cat("Ligne : Groupe ID_proteine nom_proteine\n")
  
  singletons <- find_singletons(amisX, nb_prot)
  cat("##############################################\n")
  cat("# Singletons - ", length(singletons), " proteines\n")
  for(s in 1:length(singletons)){
    cat("0 ", singletons[s], namesX[singletons[s]],"\n")
  }
  
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


getAlphabet <- function(struct_fasta){
  alpha <- list()
  for(i in 1:length(struct_fasta)) #sequences fasta
  {
    for(j in 1:length(struct_fasta[[i]])) #caracteres de la sequence fasta i
    {
      if(" " != struct_fasta[[i]][[j]]) alpha <- union(alpha, struct_fasta[[i]][[j]])
    }
  }
  return (alpha)
}

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



getCaractLocation <- function(alphabet, caract){
  for(i in 1:length(alphabet)){
    if(alphabet[i] == caract) return (i)
  }
  return (0)
}

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




ecriture_all_fichiers_fasta <- function(fullNamesAndSeq, amis){
  
  dirName <- "~/R/resultats/coupe/coupe10"
  if(!dir.exists(dirName)) dir.create(dirName)
  
  
  
  
}

