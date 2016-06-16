#Recuperation des parametres demandes dans le script
options(echo= FALSE)
args <- commandArgs(trailingOnly= TRUE)
nb_arguments <- length(args)

#variables communes aux deux methodes
cat(sprintf("Prise en compte de vos parametres"), "\n")

fic_path <- args[1]
root_dir <- dirname(dirname(fic_path))
nb_individus <- as.numeric(args[2])
min <- as.numeric(args[3])
max <- as.numeric(args[4])
hm_fic <- args[5]
data_fasta <- args[6]
#on peut vouloir les occurrences sans les fichiers fasta des groupes trouves
#il se peut donc que data_fasta == occ
occ <- args[7]
clique <- args[8]

print(args)

# nb_arguments <- 9
# fic_path <- "~/R/data/mnemo-archees_bacteries_VLD.raw"
# root_dir <- dirname(dirname(fic_path))
# print(fic_path)
# print(dirname(fic_path))
# print(dirname(dirname(fic_path)))
# nb_individus <- 2254
# min <- 5
# max <- 25
# hm_fic <- "n"
# data_fasta <- "n"
# occ <- "n"
# clique <- "n"

source(paste(root_dir, "/functions/traitement_fichier.r", sep= ""))
source(paste(root_dir, "/functions/traitementclustering.r", sep= ""))
source(paste(root_dir, "/functions/clique.r", sep= ""))

list.of.packages <- c("cluster")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/")
library(cluster)

save_file <- paste(substr(fic_path, 1, nchar(fic_path)-4), "_RET.RData", sep= "")
if(file.exists(save_file)) #donnees deja sauvegardees
{
  cat(sprintf("Recuperation des donnees precedemment chargees\n"))
  load(save_file)
  
} else #premier lancement
{
  cat(sprintf("Premier lancement du programme - chargement des donnees\n"))
  
  #CHARGEMENT DONNEES
  cat(sprintf("Chargement des donnees\n"))
  dataX <- chargement_fichier(fic_path, nb_individus)
  namesX <- getProteinNamesSYS(fic_path)

  #CHARGEMENT MATRICE DE ROBUSTESSE
  cat(sprintf("Creation de la matrice de robustesse\n"))
  matRobustesse <- build_mat_rob(dataX, nb_individus, min, max)
  
  #SAUVEGARDE DES DONNEES DANS UN FICHIER EXTERNE
  cat(sprintf("Sauvegarde de vos donnees dans un fichier externe\n"))
  save(dataX, namesX, matRobustesse, file= save_file)
}


#CREATION DU FICHIER D'OCCURRENCE
if("n" != occ)
{
  baseName <- basename(occ)
  dirName <- dirname(occ)
  data_vld <- paste(dirName, "/", substr(baseName, 1, nchar(baseName)-6), "_VLD.fasta", sep= "")
  
  #Calcul des occurrences des elements de l'alphabet VLD
  cat(sprintf("Calcul des occurrences de l'alphabet VLD\n"))
  occCaractX <- getOccCaract(data_vld)
  occCaractX <- triAlphabetDecroissant(occCaractX)
  fic_occ_name <- paste("occurrence", nb_individus, sep= "")
  ecriture_fichier_occurrence(root_dir, paste(fic_occ_name, ".txt", sep= ""), occCaractX)
  fic_occ_path <- paste(root_dir, "/resultats/occurrence/", fic_occ_name, sep= "")
  ecriture_fichier_hist_occurrence(root_dir, fic_occ_name, occCaractX)
}

#CREATION FICHIER HEATMAP DE LA MATRICE DE ROBUSTESSE
if("y" == hm_fic) #fichier heatmap demande
{
  cat(sprintf("Creation fichier heatmap de la matrice de robustesse\n"))
  file_name <- paste("matRob", nb_individus, "_", min, "-", max, ".png", sep= "")
  ecriture_fichier_heatmap(root_dir, file_name, retournementMat(matRobustesse, matRobustesse[1,1]))
}

#CREATION DE CLIQUE
if("n" != clique) #fichier clique demande
{
  cat(sprintf("Creation fichier clique de la matrice de robustesse\n"))
  
  data_clique <- paste(root_dir, "/data/clique/", sep= "")
  testDir <- paste("[ -d", data_clique,"]")
  existDir <- system(testDir)
  if(existDir >= 1) #creation du dossier s'il n'existe pas
  {
    newDir <- paste("mkdir", data_clique)
    system(newDir)
  }
  
  #creation du fichier sommets
  vertex_fic_name <- paste("sommets", nb_individus, ".txt", sep= "")
  creer_fichier_sommets(root_dir, vertex_fic_name, namesX)
  
  #creation du fichier arcs
  edge_fic_name <- paste("arcs", nb_individus, ".txt", sep= "")
  robustesse <- matRobustesse[1,1] #lien maximal entre individus
  creer_fichier_arcs(root_dir, edge_fic_name, matRobustesse, namesX, robustesse)
  
  #creation de la clique a partir de la matrice de robustesse
  if("milieu" == clique) critere <- "T"
  else if ("type" == clique) critere <- "TOP"
  clique_name <- paste("clique", nb_individus, "_", critere, ".png", sep= "")
  create_clique(root_dir, clique_name, vertex_fic_name, edge_fic_name, critere, nb_individus)
}

#RECHERCHE DES GROUPES
if(9 == nb_arguments) #Classement hierarchique
{
  cat(sprintf("Classement hierarchique\n"))
  
  methode <- args[9]
  # methode <- "ward.D2"
  
  cat(sprintf("Formation des groupes et enregistrement dans le fichier de groupes\n"))
  #HCLUST ET CUTREE
  
  typeMat <- "matRobustesse"
  clust <- hclust(as.dist(matRobustesse), method= methode, members= NULL)

  setwd(paste(root_dir, "/resultats/", sep= ""))  
  
  plot_name <- paste(typeMat, nb_individus, "_", methode, sep= "")
  png(plot_name, height= 500, width= 500)
  plot(clust, main= plot_name) #affichage du dendogramme associe
  dev.off()
  
  setwd(root_dir)
  
  cat("")
  cat("L'arbre a couper se trouve dans le dossier resultats \n")
  cat("Combien de groupes voulez-vous ?\n")
  
  nb_groupes <- scan("stdin", character(), n=1)
  if("" == nb_groupes) nb_groupes <- 1
  # nb_groupes <- 12
  
  coupe <- cutree(clust, k= nb_groupes)
  amis_coupe <- build_friend_list(coupe, nb_groupes)
  
  #ecriture des fichiers de groupes formes par cutree
  nom_fichier <- paste(typeMat, nb_individus, "_", methode, "_coupe", nb_groupes, ".txt", sep= "")
  cat(sprintf("Ecriture des fichiers de groupes\n"))
  ecriture_fichier_groupes(root_dir, "coupe", nom_fichier, amis_coupe, namesX)
  
  
  # ECRITURE FICHIERS FASTA
  if("n" != data_fasta) #fichiers fasta demandes
  {
    cat(sprintf("Creation des fichiers fasta correspondants\n"))
    data_dir <- paste(root_dir, "/data/", sep= "")
    fullNamesAndSeq_X <- getNamesAndSeq(data_dir, data_fasta)
    dir_name <- paste(root_dir, "/resultats/coupe/coupe", val_coupe, "/", sep= "")
    coupeX <- build_friend_list(coupeX, val_coupe)
    ecriture_all_fichiers_fasta(dir_name, fullNamesAndSeq_X, coupeX)
  }
  
} else #Classement non hierarchique
{
  cat(sprintf("Classement non hierarchique\n"))
  
  cat(sprintf("Formation des groupes\n"))
  robustesse_max <- length(min:max) #a faire varier eventuellement
  amisX <- get_all_friends(matRobustesse, robustesse_max)
  
  cat(sprintf("Enregistrement dans le fichier de groupes\n"))
  nb_grp <- length(amisX)
  nom_fichier <- paste("amis", nb_individus,"_", nb_grp, ".txt", sep= "")
  ecriture_fichier_groupes(root_dir, "amis", nom_fichier, amisX, namesX)
  
  # ECRITURE FICHIERS FASTA
  if("n" != data_fasta) #fichiers fasta demandes
  {
    cat(sprintf("Creation des fichiers fasta correspondants\n"))
    data_dir <- dirname(fic_path)
    fullNamesAndSeq_X <- getNamesAndSeq(data_dir, data_fasta)
    dir_name <- paste(root_dir, "/resultats/amis/amis", nb_grp, "/", sep= "")
    ecriture_all_fichiers_fasta(dir_name, fullNamesAndSeq_X, amisX)
    cat(nb_grp, sprintf(" groupes ont ete crees !\n"))
  }
}
