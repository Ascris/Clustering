setwd("~/R")
source("functions/traitement_fichier.r")
source("functions/traitementclustering.r")
source("functions/triangularite.r")
source("functions/protein_name.r")
source("functions/clique.r")
source("functions/comparaison_groupes.r")

library(cluster)
library(seqinr)

########################################
####### TEST TRIANGULATION #############
########################################

#test triangulation
# print(paste("Triplets gagnants = ",triangular(data_403, 403, 1000), "%", sep=""))
# print(paste("Triplets ultrametriques gagnants = ",triangular_ultrametric(data_403, 403, 1000), "%", sep=""))

########################################
### ECHANTILLONNAGE BEST CLUSTER #######
########################################

# nb_clusters <- get_best_clustering(data_1751, 5) #environ 4min pour tester de 2 à 5
# print(paste("Nombre de clusters : ",nb_clusters)) #meilleur trouve = 2 clusters (bonne comparaison ?)
# print(paste("Moyenne :", moyenne_echantillonnage(data_403, names403, 250, nb_clusters, 10), " pour 10 echantillons de 500 individus"))

########################################
####### DETAILS PROTEINES ALEA #########
########################################

# random_names <- names403[sample(1:403, 10)]
# print(random_names)
# print("******************************")
# get_details(random_names[2])

########################################
####### COMPTAGE MULTI CRITERE #########
########################################
# petits problemes : certaines proteines ne font pas 13-14 caracteres !
# pam1751 <- pam(data_1751, 3)

# critere <- "M"
# count_in_all_groups(pam1751, names1751, critere)

# critere <- "RG"
# count_in_all_groups(pam1751, names1751, critere)

# print(count_multi_criteres(pam1751, names1751, c("M","RG")))

########################################
############# TEST MDS #################
########################################

# pam403 <- pam(data_403, 3)

# val_k <- 3
# MDS_1 <- cmdscale(data_1751, k = 3) #environ 10sec
# titre1751 <- paste("MDS 1751 - k= ", val_k, sep="")
# plot(MDS_1, main=titre1751)

# newpam1751 <- pam(MDS_1, 3)
# ecrire_fic_clustering(newpam1751$clustering, names1751)
# print(newpam1751)
# plot(newpam1751, main="PAM from MDS 1751 - k= 3")

# ajouter_fic_heatmap(heatmap(data_403), 403)


###################################################
############## CREATION CLIQUE ####################
###################################################

# creer_fichier_sommets("proteines403.txt", names403)
# creer_fichier_arcs("robustesse403.txt", matRobustesse403, names403, 19)
# clique <- create_clique("clique403_T", "proteines403.txt","robustesse403.txt", "T", 403)

# Rprof("temps_calcul/proteines1751.out")
# creer_fichier_sommets("proteines1751.txt", names1751)
# Rprof(NULL)
# print(summaryRprof("temps_calcul/proteines1751.out"))
# 
# Rprof("temps_calcul/robustesse1751.out")
# creer_fichier_arcs("robustesse1751.txt", matRobustesse1751, names1751, 19)
# Rprof(NULL)
# print(summaryRprof("temps_calcul/robustesse1751.out"))

# Rprof("temps_calcul/clique1751.out")
# clique <- create_clique("TESTclique1751_T", "proteines1751.txt","robustesse1751.txt", "T", 1715) #environ 110 sec
# Rprof(NULL)
# print(summaryRprof("temps_calcul/clique1751.out"))



###################################################
############## Ecriture fichier groupes ###########
###################################################

# ecriture_fichier_groupes("amis403.txt", amis403, names403)
# ecriture_fichier_groupes("amis1751.txt", amis1751, names1751)


###################################################
############## Matrice de robustesse ++ ###########
###################################################

# Rprof("temps_calcul/test.out")
# super_matRob403 <- build_all_mat_rob(data_403, 403, 25, 5, 25) # environ 240 sec
# Rprof(NULL)
# print(summaryRprof("temps_calcul/test.out"))

# creer_fichier_sommets("proteines403.txt", names403)
# creer_fichier_arcs("super_robustesse403.txt", super_matRob403, names403, 21) # last param 21 = length(5:25)
# super_clique <- create_clique("proteines403.txt","super_robustesse403.txt", "TOP", 403)

# super_amis403 <- get_all_friends(super_matRob403, 21) # 1 min environ pour le 403
# print(paste("amis = ", super_amis403)) # 30 groupes

# ecriture_fichier_groupes("super_amis403.txt", super_amis403, names403)

# stabilite_groupes403 <- stabilite_groupes(data_403, 403, 1, 5, 25, names403, 5) # 12min environ pour dim= 5

# Rprof("temps_calcul/stabilite403.out")
# stabilite_groupes403 <- stabilite_groupes(data_403, 403, 25, 5, 25, names403, 10) # 48min environ
# Rprof(NULL)
# print(summaryRprof("temps_calcul/stabilite403.out"))

# Rprof("temps_calcul/stabilite1751.out")
# stabilite_groupes1751 <- stabilite_groupes(data_1751, 1751, 1, 5, 25, names1751, 1)
# Rprof(NULL)
# print(summaryRprof("temps_calcul/stabilite1751.out"))

#############################################################
############## Pistes classification plus précise ###########
#############################################################

# testHeatmap <- heatmap(data_403, distfun=as.dist, keep.dendro=T)
# print(testHeatmap$Colv)

# res = hclust(as.dist(data_403), method= "ward.D2", members = NULL)
# print(names(res))
# plot(res)


############################################################
########## Comparaison des groupes trouves  ################
############################################################


# super_matRob403 <- build_all_mat_rob(data_403, 403, 1, 5, 25) # environ 1min
# super_amis403 <- get_all_friends(super_matRob403, 10) # 1 min environ pour le 403
# 14 -> 17 groupes
# 10 -> 09 groupes
# print(paste("amis = ", super_amis403)) # 30 groupes avec 21

# solid_elements <- get_solid_elements(super_amis403, amis_coupe10)

# plop <- get_sub_lists(super_amis403)



