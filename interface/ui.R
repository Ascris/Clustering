library(shiny)
library(png)

# toutes les images integrees doivent etre dans 'www' qui doit etre avec ui.R

liste_cliques <- system("ls ../resultats/clique/ | grep '.png' ", intern= TRUE) #recuperation des cliques dans 'resultats'

shinyUI(fluidPage(
  titlePanel("Interface Shiny - Coloration de cliques"),
  
  #partie gauche de l'interface
  sidebarLayout(
    sidebarPanel(
      p("Interface concue en ",
        img(src= "bigorb.png", height= 30, width= 30),
        "avec",
        span("RStudio", style= "color:red")
      ),
      p("Apprendre a utiliser", a("shiny", href="http://shiny.rstudio.com/tutorial/")),
      br(),
      br(),
      
      fluidRow(
        column(12,
               fileInput("fichier_groupes",
                         label= "Fichier de groupes : ",
                         accept = c('text/plain'),
                         multiple= FALSE)
        )
      ),
      
      fluidRow(
        wellPanel(id= "fichier_groupes", style= "overflow-y:scroll; max-height: 600px",
                  tableOutput("liste_groupes")
        )
      )
      
    ), #fin du sidebarPanel
    #partie droite principale de l'interface
    mainPanel(
      fluidRow(
        column(6, selectInput("clique1", label= "Nom clique 1", choices= liste_cliques)),
        column(6, selectInput("clique2", label= "Nom clique 2", choices= liste_cliques))
      ),
      
      fluidRow(
        column(6, imageOutput("clique1")),
        column(6, imageOutput("clique2"))
      )
    ) #fin du mainPanel
  )# fin du sidebarLayout
))
