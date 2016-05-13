library(shiny)
library(png)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$clique1 <- renderImage({
    return(list(
      src= paste("../resultats/clique/", input$clique1, sep=""),
      contentType= "image/png",
      alt= "Clique 1",
      width= 400,
      height= 400
    ))
  }, deleteFile = FALSE)
  
  output$clique2 <- renderImage({
    return(list(
      src= paste("../resultats/clique/", input$clique2, sep=""),
      contentType= "image/png",
      alt= "Clique 2",
      width= 400,
      height= 400
    ))
  }, deleteFile = FALSE)
  
  output$liste_groupes <- renderTable({
    fic <- input$fichier_groupes
    if(!is.null(fic)) read.csv(fic$datapath)
  })
  
})
