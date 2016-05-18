setwd("~/R")
source("functions/traitement_fichier.r")
source("functions/traitementclustering.r")
source("functions/triangularite.r")
source("functions/protein_name.r")
source("functions/clique.r")
source("functions/comparaison_groupes.r")

# install.packages('seqinr')
# install.packages('VennDiagram')

library(cluster)
library(seqinr)
library(VennDiagram)

########################################
####### CHARGEMENT FICHIER 403 #########
########################################

# FICHIER 403 proteines
# data403 <- chargement_fichier("data/403_VLD_dist.raw", 403) #environ 7sec
# names403 <- getProteinNames("data/403_VLD_dist.raw")

########################################
####### CHARGEMENT FICHIER 1751 ########
########################################

# FICHIER 1751 proteines
# data_1751 <- chargement_fichier("data/1751_VLD_dist.raw", 1751) #environ 4min
# names1751 <- getProteinNames("data/1751_VLD_dist.raw")
  # tab_real_names <- get_tab_real_names("new_names1751.txt", 1751)
  # names1751 <- get_real_names(tab_real_names)


###################################################
############ ROBUSTESSE GROUPES ###################
###################################################

# matRobustesse403 <- build_mat_rob(data_403, 403, 2, 20) #environ 8 sec
# amis403 <- get_all_friends(matRobustesse403, 19) # 1 min environ pour le 403

# matRobustesse1751 <- build_mat_rob(data_1751, 1751, 2, 20) # environ 19 min pour le 1751
# amis1751 <- get_all_friends(matRobustesse1751, 19) # 120 min environ pour le 1751


###################################################
############## Matrice de robustesse ++ ###########
###################################################

# stabilite_groupes403 <- stabilite_groupes(data_403, 403, 1, 5, 25, names403, 5) # 12min environ pour dim= 5
# stabilite_groupes1751 <- stabilite_groupes(data_1751, 1751, 1, 5, 25, names1751, 5)


#######################################
############ Coupe d'arbre ############
#######################################

############# MATRICE DE DISTANCES ######################

# ##### METHODE WARD D #####
# DATAcoupe10wardD <- cutAndWrite(data_403, "DATA", "ward.D", 10, "coupe", names403)
# DATAcoupe16wardD <- cutAndWrite(data_403, "DATA", "ward.D", 16, "coupe", names403)
# ############################################################
########## Creation fichiers fasta groupes  ################
############################################################
# 
# ##### METHODE WARD D2 #####
# DATAcoupe10wardD2 <- cutAndWrite(data_403, "DATA", "ward.D2", 10, "coupe", names403)
# DATAcoupe16wardD2 <- cutAndWrite(data_403, "DATA", "ward.D2", 16, "coupe", names403)
# 
# 
# ############# MATRICE DE ROBUSTESSE INVERSEE ############
# ##### METHODE WARD D #####
# retMat <- retournementMat(matRobustesse403, matRobustesse[1,1]) #retournement de la matrice
# 
# RETcoupe10wardD <- cutAndWrite(retMat, "RET", "ward.D", 10, "coupe", names403)
# RETcoupe16wardD <- cutAndWrite(retMat, "RET", "ward.D", 16, "coupe", names403)
# 
# ##### METHODE WARD D2 #####
# RETcoupe10wardD2 <- cutAndWrite(retMat, "RET", "ward.D2", 10, "coupe", names403)
# RETcoupe16wardD2 <- cutAndWrite(retMat, "RET", "ward.D2", 16, "coupe", names403)
# 
# ############# SOUS MATRICE DE ROBUSTESSE INVERSEE ############
# subMat <- getSubAndRetMat(matRobustesse403, 403, 10) #sous matrice retournee
# 
# ##### METHODE WARD D #####
# SUBcoupe10wardD <- cutAndWrite(subMat, "SUB", "ward.D", 10, "coupe", names403)
# SUBcoupe16wardD <- cutAndWrite(subMat, "SUB", "ward.D", 16, "coupe", names403)
# 
# ##### METHODE WARD D2 #####
# SUBcoupe10wardD2 <- cutAndWrite(subMat, "SUB", "ward.D2", 10, "coupe", names403)
# SUBcoupe16wardD2 <- cutAndWrite(subMat, "SUB", "ward.D2", 16, "coupe", names403)


############################################################
########## Lecture FASTA, alphabet et occurrence ###########
############################################################

# occCaract <- getSortAndWriteOccCaract("403_seq.vld.fasta", "occCaract403.txt")

# ecriture_fichier_hist_occurrence("hist_occurrence403.png", occCaract)

############################################################
########## Creation fichiers fasta groupes  ################
############################################################

# fullNamesAndSeq_403 <- getNamesAndSeq("403_seq.fasta")
# ecriture_all_fichiers_fasta("~/R/resultats/coupe/coupe10/", fullNamesAndSeq_403, amis_coupe10)
# ecriture_all_fichiers_fasta("~/R/resultats/amis/", fullNamesAndSeq_403, amis403)

############################################################
################ Diagramme de Venn  ########################
############################################################

# On applique le diagramme de Venn sur les 3 versions de chaque groupe des coupes en 10 groupes
# du hclust cree sur les 403 proteines. On a une version en utilisant la matrice de distance vld, 
# une version sur la matrice de robustesse inversée (application de PAM X fois) et une dernière
# version sur une sous-matrice de robustesse 10% plus petite. On utilisera la methode "ward.D2"
# comme reference dans ces 3 versions de chaque groupe.

### Donnees utilisees ###
# Version matrice de distance #
# DATAcoupe10wardD2 <- cutAndWrite(data_403, "DATA", "ward.D2", 10, "coupe", names403)

#matrice de robustesse
# matRobustesse403 <- build_mat_rob(data_403, 403, 2, 20)

# Version matrice de robustesse #
# retMat <- retournementMat(matRobustesse403, matRobustesse[1,1]) #retournement de la matrice
# RETcoupe10wardD2 <- cutAndWrite(retMat, "RET", "ward.D2", 10, "coupe", names403)

# Version sous matrice de robustesse #
# subMat <- getSubAndRetMat(matRobustesse403, 403, 10) #sous matrice retournee
# SUBcoupe10wardD2 <- cutAndWrite(subMat, "SUB", "ward.D2", 10, "coupe", names403)

# DATAcoupe10wardD2 ; RETcoupe10wardD2 ; SUBcoupe10wardD2

ecriture_fichiers_venn(DATAcoupe10wardD2, RETcoupe10wardD2, SUBcoupe10wardD2, 10)












