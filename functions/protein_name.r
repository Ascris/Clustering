get_regne <- function(nom){
  return (substr(nom, nchar(nom), nchar(nom)))
}

isArchea <- function(nom){
  return ('A' == get_regne(nom))
}

isBacteria <- function(nom){
  return ('B' == get_regne(nom))
}

isEuracyote <- function(nom){
  return ('E' == get_regne(nom))
}

#nom raccourci de la proteine
get_genre <- function(nom){
  return (substr(nom, 1, 6))
}

#reverse gyrase ou topoisomerase
get_type <- function(nom){
  if(13 == nchar(nom)) return (substr(nom, 7, 8))
  return (substr(nom, 7, 9))
}

isTopoisomerase <- function(nom){
  return("TOP" == get_type(nom))
}

isReverseGyrase <- function(nom){
  return("RG" == get_type(nom))
}

get_sousfamille <- function(nom){
  return(substr(nom, nchar(nom)-4, nchar(nom)-4))
}

# milieu : thermophile, mezophile, psychrophile, halophile
# Thermophile = aime les temperatures elevees (+80°C)
# Mezophile = aime les temperatures moderees (20-40°C)
# Psychrophile = aime les temperatures basses (-15°C)
# Halophile = aime les milieux sales
get_milieu <- function(nom){
  return (substr(nom, nchar(nom)-2, nchar(nom)-2))
}

isThermophile <- function(nom){
  return ("T" == get_milieu(nom))
}

isMezophile <- function(nom){
  return ("M" == get_milieu(nom))
}

isPsychrophile <- function(nom){
  return ("P" == get_milieu(nom))
}

isHalophile <- function(nom){
  return ("H" == get_milieu(nom))
}

get_full_regne <- function(regne){
  return (switch(regne,"A"="Archea","B"="Bacteria","E"="Eucaryote"))
}

get_full_milieu <- function(milieu){
  return (switch(milieu,"T"="Thermophile","M"="Mezophile","P"="Psychrophile", "H"="Halophile", "X"="???"))
}

get_full_type <- function(type){
  return (switch(type,"TOP"="Topoisomerase","RG"="Reverse Gyrase"))
}

#retourne les informations relatives a la proteine "nom"
get_details <- function(nom){
  print(paste(nom, " possede les caracteristiques suivantes :"))
  print(paste("regne = ", get_full_regne(get_regne(nom))))
  print(paste("Milieu = ", get_full_milieu(get_milieu(nom))))
  print(paste("type = ", get_full_type(get_type(nom))))
  print(paste("Sous famille = ", get_sousfamille(nom)))
  print(paste("Genre = ", get_genre(nom)))
  
}

#which_critere retourne une chaine correspondant a la nature du critere (regne, type, ou milieu)
which_critere <- function(critere){
  x <- switch(critere, 
              "A"="regne", "B"="regne", "E"="regne",
              "TOP"="type", "RG"="type",
              "T"="milieu","M"="milieu","P"="milieu","H"="milieu")
  
  return (x)
}

# Retourne 1 si les elements correspondent, 0 sinon
testRegne <- function(prot, critere){
  if(critere == get_regne(prot)) return (1)
  else return (0)
}

# Retourne 1 si les elements correspondent, 0 sinon
testType <- function(prot, critere){
  if(critere == get_type(prot)) return (1)
  else return (0)
}

# Retourne 1 si les elements correspondent, 0 sinon
testMilieu <- function(prot, critere){
  if(critere == get_milieu(prot)) return (1)
  else return (0)
}

#retourne le nom complet du critere
get_full_critere <- function(critere){
  full <- switch(which_critere(critere),
                 "regne"=get_full_regne(critere),
                 "milieu"=get_full_milieu(critere),
                 "type"=get_full_type(critere))
  return (full)
}