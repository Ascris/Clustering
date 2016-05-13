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
####### CHARGEMENT FICHIER 403 #########
########################################

# FICHIER 403 proteines
# data403 <- chargement_fichier("data/403_VLD_dist.raw", 403) #environ 7sec
# data_403 <- remplacementNApar0(data_403)
# names403 <- getProteinNames("data/403_VLD_dist.raw")

########################################
####### CHARGEMENT FICHIER 1751 ########
########################################

# FICHIER 1751 proteines
# data_1751 <- chargement_fichier("data/1751_VLD_dist.raw", 1751) #environ 4min
# data_1751 <- remplacementNApar0(data_1751)
# names1751 <- getProteinNames("data/1751_VLD_dist.raw")
  # tab_real_names <- get_tab_real_names("new_names1751.txt", 1751)
  # names1751 <- get_real_names(tab_real_names)


########################################
####### TEST TRIANGULATION #############
########################################

#test triangulation
# print(paste("Triplets gagnants = ",triangular(data_403, 403, 1000), "%", sep=""))
# print(paste("Triplets ultrametriques gagnants = ",triangular_ultrametric(data_403, 403, 1000), "%", sep=""))


########################################
### ECHANTILLONNAGE BEST CLUSTER #######
########################################

# nb_clusters <- get_best_clustering(data_1751, 2) #environ 4min pour tester de 2 à 5
# print(paste("Nombre de clusters : ",nb_clusters)) #meilleur trouve = 2 clusters (bonne comparaison ?)
# print(paste("Moyenne :", moyenne_echantillonage(data_403, names403, 250, nb_clusters, 10), " pour 10 echantillons de 500 individus"))

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
############ ROBUSTESSE GROUPES ###################
###################################################

# matRobustesse403 <- build_mat_rob(data_403, 403, 2, 20) #environ 8 sec
# print(matRobustesse403[1:20,1:20])
# amis403 <- get_all_friends(matRobustesse403, 19) # 1 min environ pour le 403
# print(paste("amis = ", amis403))

# matRobustesse1751 <- build_mat_rob(data_1751, 1751, 2, 20) # environ 19 min pour le 1751
# print(matRobustesse1751[1:20,1:20])
# amis1751 <- get_all_friends(matRobustesse1751, 19) # 120 min environ pour le 1751
# print(paste("amis = ",amis1751))

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

#######################################
############ Coupe d'arbre ############
#######################################

############# MATRICE DE DISTANCES ######################

##### METHODE WARD D #####
# retMat <- data_403
# clust_retMat <- hclust(as.dist(retMat), method = "ward.D", members= NULL)
# plot(clust_retMat)
# 
# coupe_16 <- cutree(clust_retMat, k= 16)
# amis_coupe16 <- build_friend_list(coupe_16, 16)
# ecriture_fichier_groupes("coupe", "DATAmat403ret19_wardD_coupe16.txt", amis_coupe16, names403)

# coupe_10 <- cutree(clust_retMat, k= 10)
# amis_coupe10 <- build_friend_list(coupe_10, 10)
# ecriture_fichier_groupes("coupe", "DATAmat403ret19_wardD_coupe10.txt", amis_coupe10, names403)
# 
# coupe_3 <- cutree(clust_retMat, k= 3)
# amis_coupe3 <- build_friend_list(coupe_3, 3)
# ecriture_fichier_groupes("coupe", "DATAmat403ret19_wardD_coupe3.txt", amis_coupe3, names403)
# 
# coupe_4 <- cutree(clust_retMat, k= 4)
# amis_coupe4 <- build_friend_list(coupe_4, 4)
# ecriture_fichier_groupes("coupe", "DATAmat403ret19_wardD_coupe4.txt", amis_coupe4, names403)

##### METHODE WARD D2 #####
# retMat <- data_403
# clust_retMat <- hclust(as.dist(retMat), method = "ward.D2", members= NULL)
# plot(clust_retMat)
# 
# coupe_16 <- cutree(clust_retMat, k= 16)
# amis_coupe16 <- build_friend_list(coupe_16, 16)
# ecriture_fichier_groupes("coupe", "DATAmat403ret19_wardD2_coupe16.txt", amis_coupe16, names403)

# 
# coupe_10 <- cutree(clust_retMat, k= 10)
# amis_coupe10 <- build_friend_list(coupe_10, 10)
# ecriture_fichier_groupes("coupe", "DATAmat403ret19_wardD2_coupe10.txt", amis_coupe10, names403)
# 
# coupe_3 <- cutree(clust_retMat, k= 3)
# amis_coupe3 <- build_friend_list(coupe_3, 3)
# ecriture_fichier_groupes("coupe", "DATAmat403ret19_wardD2_coupe3.txt", amis_coupe3, names403)
# 
# coupe_4 <- cutree(clust_retMat, k= 4)
# amis_coupe4 <- build_friend_list(coupe_4, 4)
# ecriture_fichier_groupes("coupe", "DATAmat403ret19_wardD2_coupe4.txt", amis_coupe4, names403)


############# MATRICE DE ROBUSTESSE INVERSEE ############
##### METHODE WARD D #####
# retMat <- retournementMat(matRobustesse403, matRobustesse[1,1])
# clust_retMat <- hclust(as.dist(retMat), method = "ward.D", members= NULL)
# plot(clust_retMat)
# 
# coupe_16 <- cutree(clust_retMat, k= 16)
# amis_coupe16 <- build_friend_list(coupe_16, 16)
# ecriture_fichier_groupes("coupe", "mat403ret19_wardD_coupe16.txt", amis_coupe16, names403)

# 
# coupe_10 <- cutree(clust_retMat, k= 10)
# amis_coupe10 <- build_friend_list(coupe_10, 10)
# ecriture_fichier_groupes("coupe", "mat403ret19_wardD_coupe10.txt", amis_coupe10, names403)
# 
# coupe_3 <- cutree(clust_retMat, k= 3)
# amis_coupe3 <- build_friend_list(coupe_3, 3)
# ecriture_fichier_groupes("coupe", "mat403ret19_wardD_coupe3.txt", amis_coupe3, names403)
# 
# coupe_4 <- cutree(clust_retMat, k= 4)
# amis_coupe4 <- build_friend_list(coupe_4, 4)
# ecriture_fichier_groupes("coupe", "mat403ret19_wardD_coupe4.txt", amis_coupe4, names403)

##### METHODE WARD D2 #####
# retMat <- retournementMat(matRobustesse403, matRobustesse[1,1])
# clust_retMat <- hclust(as.dist(retMat), method = "ward.D2", members= NULL)
# plot(clust_retMat)
# 
# coupe_16 <- cutree(clust_retMat, k= 16)
# amis_coupe16 <- build_friend_list(coupe_16, 16)
# ecriture_fichier_groupes("coupe", "mat403ret19_wardD2_coupe16.txt", amis_coupe16, names403)

# 
# coupe_10 <- cutree(clust_retMat, k= 10)
# amis_coupe10 <- build_friend_list(coupe_10, 10)
# ecriture_fichier_groupes("coupe", "mat403ret19_wardD2_coupe10.txt", amis_coupe10, names403)
# 
# coupe_3 <- cutree(clust_retMat, k= 3)
# amis_coupe3 <- build_friend_list(coupe_3, 3)
# ecriture_fichier_groupes("coupe", "mat403ret19_wardD2_coupe3.txt", amis_coupe3, names403)
# 
# coupe_4 <- cutree(clust_retMat, k= 4)
# amis_coupe4 <- build_friend_list(coupe_4, 4)
# ecriture_fichier_groupes("coupe", "mat403ret19_wardD2_coupe4.txt", amis_coupe4, names403)

#####################################################
## Extraction individus et tests stabilité groupes ##
#####################################################

# nb_prot <- 403
# percent <- 10 # pourcentage a enlever
# remaining_ones <- sort(sample(1:nb_prot, nb_prot - (nb_prot*percent/100)))
# print(remaining_ones)
# print(paste("taille de l'echantillon : ", length(remaining_ones)))


# Methode ward.D
# retMat <- retournementMat(matRobustesse403[remaining_ones, remaining_ones], matRobustesse[1,1])
# clust_retMat <- hclust(as.dist(retMat), method = "ward.D", members= NULL)
# plot(clust_retMat)
# 
# coupe_16 <- cutree(clust_retMat, k= 16)
# amis_coupe16 <- build_friend_list(coupe_16, 16)
# ecriture_fichier_groupes("coupe", "SUBmat403ret19_wardD_coupe16.txt", amis_coupe16, names403)
# 
# coupe_10 <- cutree(clust_retMat, k= 10)
# amis_coupe10 <- build_friend_list(coupe_10, 10)
# ecriture_fichier_groupes("coupe", "SUBmat403ret19_wardD_coupe10.txt", amis_coupe10, names403)
# 
# Methode ward.D2
# retMat <- retournementMat(matRobustesse403[remaining_ones, remaining_ones], matRobustesse[1,1])
# clust_retMat <- hclust(as.dist(retMat), method= "ward.D2", members= NULL)
# plot(clust_retMat)
# 
# coupe_16 <- cutree(clust_retMat, k= 16)
# amis_coupe16 <- build_friend_list(coupe_16, 16)
# ecriture_fichier_groupes("coupe", "SUBmat403ret19_wardD2_coupe16.txt", amis_coupe16, names403)
# 
# coupe_10 <- cutree(clust_retMat, k= 10)
# amis_coupe10 <- build_friend_list(coupe_10, 10)
# ecriture_fichier_groupes("coupe", "SUBmat403ret19_wardD2_coupe10.txt", amis_coupe10, names403)

############################################################
########## Lecture FASTA, alphabet et occurrence ###########
############################################################

# occCaract <- getOccCaract("403_seq.vld.fasta") #environ 2h
# new_occCaract <- triAlphabet(occCaract)
# ecriture_fichier_occurrence("occCaract403.txt", new_occCaract)


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


############################################################
########## Creation fichiers fasta groupes  ################
############################################################


fullNamesAndSeq_403 <- getNamesAndSeq("403_seq.fasta")

create_fichier_fasta("", "", fullNamesAndSeq_403, "")








