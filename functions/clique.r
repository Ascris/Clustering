current_dir <- getwd()
source(paste(current_dir, "/functions/traitement_fichier.r", sep= ""))
source(paste(current_dir, "/functions/protein_name.r", sep= ""))

list.of.packages <- c("igraph")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/")

library(igraph)


#' Cree et retourne la clique formee des proteines du 'fic_sommets' liees par le 'fic_arcs
#' la clique retournee est coloree en fonction du 'critere' (regne, type ou milieu) et correspond
#' au regne (A, B ou E), a le type (TOP ou RG) ou au milieu (T, M, P, H)
#'
#' @param root_dir : dossier racine de l'application (chaine de caracteres) 
#' @param clique_name : nom du fichier clique desire (chaine de caracteres)
#' @param fic_sommets : fichier contenant les proteines (ex : "AcamarTOPa_M_B" "B" "TOP" "M") (chaine de caracteres)
#' @param fic_arcs : fichier contenant la robustesse entre les proteines (ex : "AcibooRGz_T_A" "AcamarTOPa_M_B" 20) (chaine de caracteres)
#' @param critere : selection de la coloration du graphe ("regne", "type" ou "milieu") (chaine de caracteres)
#' @param nb_prot : nombre de proteines prises en compte (entier)
#'
#' @return plot de la clique avec sa legende
#' @export Cree une clique et sa legende et l'enregistre en .png dans 'resultats/clique/'
#'
#' @examples clique <- create_clique("~/R/data/", "clique403_T", "proteines403.txt","robustesse403.txt", "T", 403)
create_clique <- function(root_dir, clique_name, fic_sommets, fic_arcs, critere, nb_prot){
  currentDir <- getwd()
  data <- paste(root_dir, "/data/clique/", sep= "")
  setwd(data)
  
  prot <- read.table(fic_sommets, sep=" ", header= TRUE, row.names= NULL)
  robustesse <- read.table(fic_arcs, sep=" ", header= TRUE, stringsAsFactors= FALSE, row.names= NULL)
  
  clique <- graph.data.frame(robustesse, directed= FALSE, vertices= prot)
  
  x <- switch(critere,
              "A"="regne", "B"="regne", "E"="regne",
              "TOP"="type", "RG"="type",
              "T"="milieu","M"="milieu","P"="milieu","H"="milieu")

  two_elements <- TRUE
  
  if("milieu" == x) #milieu = T, M, P ou H
  {
    print("dans cas : milieu")
    legende_couleur <- c("T","M","P","H","?")
    couleur_prot <- get.vertex.attribute(clique, "type")
    couleur_prot[is.na(couleur_prot)] <- "?"
    couleurs <- c("Red", "Green", "Blue", "Yellow", "Pink")
    for(i in 1:length(couleur_prot)){
      if("T" == couleur_prot[i]) couleur_prot[i] <- couleurs[1]
      else if("M" == couleur_prot[i] || "X" == couleur_prot[i]) couleur_prot[i] <- couleurs[2]
      else if("P" == couleur_prot[i]) couleur_prot[i] <- couleurs[3]
      else if("H" == couleur_prot[i]) couleur_prot[i] <- couleurs[4]
      else couleur_prot[i] <- couleurs[5]
    }
    
    legende_taille <- c("A", "B", "E", "?")
    taille_points <- get.vertex.attribute(clique, "nom_prot")
    taille_points[is.na(taille_points)] <- "?"
    tailles <- c(15, 10, 8, 3)
    for(i in 1:length(taille_points)){
      if("A" == taille_points[i]) taille_points[i] <- tailles[1]
      else if("B" == taille_points[i]) taille_points[i] <- tailles[2]
      else if("E" == taille_points[i]) taille_points[i] <- tailles[3]
      else taille_points[i] <- tailles[4]
    }
    titre_legende_1 <- "Milieu"
    titre_legende_2 <- "Regne"
    
  } else
    {
      if("type" == x) #type = TOP ou RG
      {
        print("dans cas : type")
        legende_couleur <- c("TOP", "RG", "?")
        couleur_prot <- get.vertex.attribute(clique, "regne")
        couleur_prot[is.na(couleur_prot)] <- "?"
        couleurs <- c("Red", "Green", "Blue")
        for(i in 1:length(couleur_prot)){
          if("TOP" == couleur_prot[i]) couleur_prot[i] <- couleurs[1]
          else if("RG" == couleur_prot[i]) couleur_prot[i] <- couleurs[2]
          else couleur_prot[i] <- couleurs[3]
        }
        
        legende_taille <- c("A", "B", "E", "?")
        taille_points <- get.vertex.attribute(clique, "nom_prot")
        taille_points[is.na(taille_points)] <- "?"
        tailles <- c(15, 10, 8, 3)
        for(i in 1:length(taille_points)){
          if("A" == taille_points[i]) taille_points[i] <- tailles[1]
          else if("B" == taille_points[i]) taille_points[i] <- tailles[2]
          else if("E" == taille_points[i]) taille_points[i] <- tailles[3]
          else taille_points[i] <- tailles[4]
        }
        titre_legende_1 <- "Type"
        titre_legende_2 <- "Regne"
        
      } else # "regne" == x par defaut (A, B ou E)
        {
          print("dans cas : regne ou autre")
          two_elements <- FALSE
          legende_couleur <- c("A", "B", "E", "?")
          couleur_prot <- get.vertex.attribute(clique, "nom_prot")
          couleur_prot[is.na(couleur_prot)] <- "?"
          couleurs <- c("Red", "Green", "Blue", "Yellow")
          for(i in 1:length(couleur_prot)){
            if("A" == couleur_prot[i]) couleur_prot[i] <- couleurs[1]
            else if("B" == couleur_prot[i]) couleur_prot[i] <- couleurs[2]
            else if("E" == couleur_prot[i]) couleur_prot[i] <- couleurs[3]
            else couleur_prot[i] <- couleurs[4]
          }
        }
        titre_legende_1 <- "Regne"
    
  }
  
  taille_arcs <- get.edge.attribute(clique, "force")
  
  #CREATION FICHIER DESSIN
  res_clique <- paste(root_dir, "/resultats/clique/", sep= "")
  testDir <- paste("[ -d", res_clique,"]")
  existDir <- system(testDir)
  if(existDir >= 1) #creation du dossier s'il n'existe pas
  {
    newDir <- paste("mkdir", res_clique)
    system(newDir)
  }
  setwd(res_clique)
  
  fic_name <- clique_name
  new_image <- get_next_filename(fic_name) # nom de fichier non existant
  png(filename= new_image, width=500, height=500)
  
  layout(matrix(c(1,2,3), 3, byrow = TRUE), height= c(7,1,2)) # positionnement des elements de l'image
  
  if(two_elements)
  {
    plot_clique <- plot(clique,
                        vertex.color= couleur_prot,
                        edge.width= taille_arcs,
                        edge.color= "black",
                        vertex.label= NA,
                        vertex.size= as.matrix(as.numeric(taille_points)))
    
  } else {
    plot_clique <- plot(clique,
                        vertex.color= couleur_prot,
                        edge.width= taille_arcs,
                        edge.color= "black",
                        vertex.label= NA)
  }
  par(mar=c(0,0,0,0))
  plot.new()
  legend("center", "groups", legend= legende_couleur, pch= 16, col= couleurs, ncol=5, title= titre_legende_1)
  
  if(two_elements){
    plot.new()
    legend("center", "groups", legend= legende_taille, pch= 16, cex = 2, pt.cex= tailles/3, ncol= 4, title= titre_legende_2)
  }
    
  dev.off()
  #FIN CREATION FICHIER DESSIN
  
  setwd(currentDir)
  return (plot_clique)
}
