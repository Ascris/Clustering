root_dir <- getwd()

source(paste(root_dir, "/functions/traitement_fichier.r", sep= ""))
source(paste(root_dir, "/functions/traitementclustering.r", sep= ""))

library(cluster)

#Recuperation des parametres demandes dans le script
options(echo= FALSE)
args <- commandArgs(trailingOnly= TRUE)
nb_arguments <- length(args)

#variables communes aux deux methodes
cat(sprintf("Prise en compte de vos parametres"), "\n")
fic_path <- args[1]
nb_individus <- as.numeric(args[2])
min <- as.numeric(args[3])
max <- as.numeric(args[4])
hm_fic <- args[5]
data_fasta <- args[6]

#CHARGEMENT DONNEES
cat(sprintf("Chargement des donnees"), "\n")
dataX <- chargement_fichier(fic_path, nb_individus)
namesX <- getProteinNames(fic_path)

#CHARGEMENT MATRICE DE ROBUSTESSE
cat(sprintf("Creation de la matrice de robustesse"), "\n")
matRobustesse <- build_mat_rob(dataX, nb_individus, min, max)

if("y" == hm_fic) #fichier heatmap demande
{
  cat(sprintf("Creation fichier heatmap de la matrice de robustesse"), "\n")
  file_name <- paste("matRob", nb_individus, "_", min, "-", max, ".png", sep= "")
  ecriture_fichier_heatmap(root_dir, file_name, matRobustesse)
}

if(8 == nb_arguments) #Classement hierarchique
{
  cat(sprintf("Classement hierarchique"), "\n")
  
  methode <- args[7]
  val_coupe <- as.numeric(args[8])
  
  cat(sprintf("Formation des groupes et enregistrement dans le fichier de groupes"), "\n")
  #HCLUST ET CUTREE
  coupeX <- cutAndWrite(dataX, "Coupe", methode, val_coupe, root_dir, "coupe", namesX)
  
  # ECRITURE FICHIERS FASTA
  if("n" != data_fasta) #fichiers fasta demandes
  {
    cat(sprintf("Creation des fichiers fasta correspondants"), "\n")
    data_dir <- paste(root_dir, "/data/", sep= "")
    fullNamesAndSeq_X <- getNamesAndSeq(data_dir, data_fasta)
    dir_name <- paste(root_dir, "/resultats/coupe/coupe", val_coupe, "/", sep= "")
    coupeX <- build_friend_list(coupeX, val_coupe)
    ecriture_all_fichiers_fasta(dir_name, fullNamesAndSeq_X, coupeX)
  }
  
} else #Classement non hierarchique
{
  cat(sprintf("Classement non hierarchique"), "\n")
  
  cat(sprintf("Formation des groupes"), "\n")
  robustesse_max <- length(min:max)
  amisX <- get_all_friends(matRobustesse, robustesse_max)
  
  cat(sprintf("Enregistrement dans le fichier de groupes"), "\n")
  nb_prot <- length(namesX)
  nb_grp <- length(amisX)
  nom_fichier <- paste("amis", nb_prot,"_", nb_grp, ".txt", sep= "")
  ecriture_fichier_groupes(root_dir, "amis", nom_fichier, amisX, namesX)
  
  # ECRITURE FICHIERS FASTA
  if("n" != data_fasta) #fichiers fasta demandes
  {
    cat(sprintf("Creation des fichiers fasta correspondants"), "\n")
    data_dir <- paste(root_dir, "/data/", sep= "")
    fullNamesAndSeq_X <- getNamesAndSeq(data_dir, data_fasta)
    dir_name <- paste(root_dir, "/resultats/amis/amis", length(amisX), "/", sep= "")
    ecriture_all_fichiers_fasta(dir_name, fullNamesAndSeq_X, amisX)
  }
}

cat(sprintf("Fin du programme"), "\n")





