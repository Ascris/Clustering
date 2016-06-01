#!/bin/bash

ROOT=$(pwd)
DATA_DIR="$ROOT/data/"
SCRIPT_DIR="$ROOT/scripts/"

echo "Bonjour et bienvenue dans ce programme de clusterisation !"
echo "Assurez-vous d'avoir bien mis votre matrice de distances (fichier .raw) dans le dossier 'data'."
echo ""

read -p "Veuillez indiquer le nom de votre matrice de distances (fichier .raw) : " FIC_NAME

FIC_PATH=$DATA_DIR$FIC_NAME
NB_ESSAIS=1

while [ ! -f $FIC_PATH ] && [ $NB_ESSAIS -lt 5 ]
do
  let NB_ESSAIS=NB_ESSAIS+1
  echo "Desole, fichier inexistant dans $DATA_DIR"
  echo ""

  read -p "Veuillez indiquer le nom de votre matrice de distances (fichier .raw) : " FIC_NAME
  
  echo "Vous avez saisi : $FIC_NAME"
  FIC_PATH=$DATA_DIR$FIC_NAME
  echo ""
done

if [ $NB_ESSAIS -eq 10 ] ; then
  echo "Veuillez placer votre fichier dans le dossier data avant de relancer le programme. A BIENTOT !"
  exit
else
  echo "Fichier de donnees existant dans $DATA_DIR !"
  echo ""
fi

#A ce stade, le fichier a ete trouve

#comptage des proteines du fichier lu
NB_INDIVIDUS=$(wc -l $FIC_PATH | cut -f1 -d' ')
echo "$NB_INDIVIDUS individus ont ete reconnus dans le fichier !"

read -p "Combien de groupes PAM minimum ? " MIN
read -p "Combien de groupes PAM maximum ? " MAX

MAIN_R=$SCRIPT_DIR"/main_script.R"

read -p "Souhaitez-vous une heatmap de la matrice de robustesse ? (y/n) " HEATMAP

#Soit on ne veut pas les fichiers et FASTA="n", soit FASTA prend le nom du fichier de donnees
read -p "Souhaitez-vous les fichiers fasta des groupes trouves en sortie ? (y/n) " FASTA

echo "Les approches proposees dans ce programme sont :" ; echo ""
echo "  - Hierarchique avec 'hclust' (0)" ; echo ""
echo "  - Non hierarchique avec reseau d'amis (1)"
read -p "Quel module desirez-vous (0 ou 1) ? " MODULE ; echo ""

if [ $MODULE -eq 1 ]
then
  echo "Non hierarchique selectionne" ; echo ""
  
  if [ "y" == $FASTA ]
  then
    read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " FASTA
  
    while [ "n" != $FASTA ] && [ ! -e $DATA_DIR$FASTA ]
    do
      echo "Fichier non trouve dans $DATA_DIR (taper \"n\" pour ne plus avoir les fichiers fasta)"
      read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " FASTA
    done
    
    if [ "n" != $FASTA ]
    then
      FASTA=$DATA_DIR$FASTA
    fi
  fi
  
  read -p "Voulez-vous les occurrences des caracteres de l'alphabet VLD ? (y/n)" OCC
  if [ "n" != $OCC ]
  then
    if [ "n" == $FASTA ]
    then
      read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " OCC
      
      while [ "n" != $OCC ] && [ ! -e $DATA_DIR$OCC ]
      do
        echo "Fichier non trouve dans $DATA_DIR (taper \"n\" pour ne plus avoir les fichiers fasta)"
        read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " OCC
      done
    
      if [ "n" != $OCC ]
      then
       OCC=$DATA_DIR$OCC
      fi
    else
      OCC=$FASTA
    fi
  fi
  
  echo "Lancement du programme avec vos parametres."
  echo "Veuillez patienter..."
  Rscript $MAIN_R $FIC_PATH $NB_INDIVIDUS $MIN $MAX $HEATMAP $FASTA $OCC #transfert des choix de l'utilisateur au script R
  
else
  echo "Hierarchique selectionne" ; echo ""
  
  NB_ESSAIS_MET=0
  METHODE=""
  while [ "ward.D" != "$METHODE" ] && [ "ward.D2" != "$METHODE" ] && [ "complete" != "$METHODE" ] &&
        [ "single" != "$METHODE" ] && [ "average" != "$METHODE" ] && [ "mcquitty" != "$METHODE" ] &&
        [ "median" != "$METHODE" ] && [ "centroid" != "$METHODE" ] && [ $NB_ESSAIS_MET -ne 5 ]
  do
    echo "Vous avez le choix parmi ces methodes : "
      echo "-\"ward.D\" "
      echo "-\"ward.D2\""
      echo "-\"complete\""
      echo "-\"single\""
      echo "-\"average (UPGMA)\""
      echo "-\"mcquitty\""
      echo "-\"median\""
      echo "-\"centroid\""
    read -p "Quelle methode utiliser ? " METHODE ; echo ""
    let NB_ESSAIS_MET=NB_ESSAIS_MET+1
  done
  
  if [ $NB_ESSAIS_MET -eq 5 ]
  then
    echo "Methode non reconnue, \"ward.D2\" choisi par defaut"
    METHODE="ward.D2"
  fi
  
  read -p "Combien de groupes voulez-vous ? " NB_GROUPES ; echo ""
  
  if [ "y" == $FASTA ]
  then
    read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " FASTA
  
    while [ "n" != $FASTA ] && [ ! -e $DATA_DIR$FASTA ]
    do
      echo "Fichier non trouve dans $DATA_DIR (taper \"n\" pour ne plus avoir les fichiers fasta)"
      read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " FASTA
    done
    
    if [ "n" != $FASTA ]
    then
      FASTA=$DATA_DIR$FASTA
    fi
  fi
  
  read -p "Voulez-vous les occurrences des caracteres de l'alphabet VLD ? (y/n)" OCC
  if [ "n" != $OCC ]
  then
    if [ "n" == $FASTA ]
    then
      read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " OCC
      
      while [ "n" != $OCC ] && [ ! -e $DATA_DIR$OCC ]
      do
        echo "Fichier non trouve dans $DATA_DIR (taper \"n\" pour ne plus avoir les fichiers fasta)"
        read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " OCC
      done
    
      if [ "n" != $OCC ]
      then
       OCC=$DATA_DIR$OCC
      fi
    else
      OCC=$FASTA
    fi
  fi
  
  echo "Lancement du programme avec vos parametres."
  echo "Veuillez patienter..."
  Rscript $MAIN_R $FIC_PATH $NB_INDIVIDUS $MIN $MAX $HEATMAP $FASTA $OCC $METHODE $NB_GROUPES #transfert des choix de l'utilisateur au script R
fi

echo "Fin du programme."