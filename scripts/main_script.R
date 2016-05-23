source("functions/traitement_fichier.r")
source("functions/traitementclustering.r")

library(cluster)

# Etapes :
#   Chargement fichier .raw
#       --> chargement donnees
#       --> creation matrice robustesse
#         --> groupes PAM pour matRob de min à max
#
#     *** Demande si heatmap sur matRob ***
#       --> Heatmap sur matRob
#
#   Methode 1 : Hclust et cutree
#     *** demande parametres : groupes voulus, methode a utiliser, nom_fic
#       
#
#   Methode 2 : Recherche d'amis
#     Avec la matRob et length(min:max) -> get_all_friends()
#     
#   --> creation fichiers groupes
#     *** fichiers fasta à demander ***

options(echo= TRUE)
args <- commandArgs(trailingOnly= TRUE)
nb_arguments <- length(args)
print(args)
print(paste("Nombre d'arguments : ", nb_arguments))

#variables communes aux deux methodes
fic_path <- args[1] ; 
root_dir <- dirname(dirname((fic_path)))
nb_individus <- as.numeric(args[2])
min <- as.numeric(args[3])
max <- as.numeric(args[4])
hm_fic <- args[5]
fasta_fic <- args[6]

#CHARGEMENT DONNEES
dataX <- chargement_fichier(fic_path, nb_individus)
namesX <- getProteinNames(fic_path)


#CHARGEMENT MATRICE DE ROBUSTESSE
matRobustesse <- build_mat_rob(dataX, nb_individus, min, max)


if(8 == nb_arguments) #Clustering hierarchique
{
  print("Clustering hierarchique")
  
  methode <- args[7]
  val_coupe <- as.numeric(args[8])
  
  #HCLUST ET CUTREE
  coupeX <- cutAndWrite(dataX, "COUPE", methode, val_coupe, root_dir, "coupe", namesX)
  
  # ECRITURE FICHIERS FASTA
  if("y" == fasta_fic)
  {
    fullNamesAndSeq_X <- getNamesAndSeq("403_seq.fasta") ##############NOM VARIABLE A CHANGER###############
    dir_name <- paste(root_dir, "resultats/coupe/coupe", val_coupe, "/", sep= "")
    ecriture_all_fichiers_fasta(dir_name, fullNamesAndSeq_X, coupeX)
  }
  
} else #Clustering non hierarchique
{
  print("Clustering non hierarchique")
  
  robustesse_max <- length(min:max)
  amisX <- get_all_friends(matRobustesse, robustesse_max)
  
  # ECRITURE FICHIERS FASTA
  if("y" == fasta_fic)
  {
    fullNamesAndSeq_X <- getNamesAndSeq("403_seq.fasta") ##############NOM VARIABLE A CHANGER###############
    dir_name <- paste(root_dir, "resultats/amis/amis", length(amisX), "/", sep= "")
    ecriture_all_fichiers_fasta(dir_name, fullNamesAndSeq_X, coupeX)
  }
}







