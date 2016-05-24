#!/bin/bash

# Etapes :
#   Chargement fichier .raw
#     *** Demande nom fichier ***
#       --> chargement donnees
#       --> creation matrice robustesse
#
#     *** Demande si heatmap sur matRob ***
#       --> Heatmap sur matRob
#
#   Hclust et cutree
#     *** demande parametres : groupes voulus, methode a utiliser, nom_fic
#       --> creation fichiers groupes
#       *** fichiers fasta à demander ***

CURRENT_DIR=$(pwd)
cd .. #repertoire racine du programme
ROOT=$(pwd)
DATA_DIR="$ROOT/data/"

echo "Bonjour et bienvenue dans ce programme de clusterisation !"
echo "Ce projet a ete realise dans le cadre dun stage en M1 Informatique"

echo "ATTENTION : Assurez-vous d'avoir bien mis votre matrice de distances (fichier .raw) dans le dossier 'data'."
echo ""

##TMP### read -p "Veuillez indiquer le nom de votre matrice de distances (fichier .raw) : " FIC_NAME
FIC_NAME="403_VLD_dist.raw"

FIC_PATH=$DATA_DIR$FIC_NAME
NB_ESSAIS=1

while [ ! -f $FIC_PATH ] && [ $NB_ESSAIS -lt 10 ]
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

###TMP### read -p "Combien de groupes PAM minimum ? " MIN
###TMP### read -p "Combien de groupes PAM maximum ? " MAX
MIN=5
MAX=6

MAIN_R=$CURRENT_DIR"/main_script.R"

###TMP### read -p "Souhaitez-vous une heatmap de la matrice de robustesse ? (y/n) " HEATMAP

#Soit on ne veut pas les fichiers et FASTA="n", soit FASTA prend le nom du fichier de donnees
###TMP### read -p "Souhaitez-vous les fichiers fasta des groupes trouves en sortie ? (y/n) " FASTA
HEATMAP="y"
FASTA="y"

echo "Les approches proposees dans ce programme sont :" ; echo ""
echo "  - Hierarchique avec 'hclust' (0)" ; echo ""
echo "  - Non hierarchique avec reseau d'amis (1)"
###TMP### read -p "Quel module desirez-vous (0 ou 1) ? " MODULE ; echo ""

MODULE=0

if [ $MODULE -eq 1 ]
then
  echo "Non hierarchique selectionne" ; echo ""
  
  if [ "y" == $FASTA ]
  then
    # read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " DATA_FASTA
    FASTA="403_seq.fasta"
  
    while [ "n" != $FASTA ] && [ ! -f $FASTA_PATH ]
    do
      echo "Fichier non trouve dans $DATA_DIR (taper \"n\" pour ne pas avoir les fichiers fasta)"
      # read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " DATA_FASTA
      FASTA="403_seq.fasta"
    done
  fi

  Rscript $MAIN_R $FIC_PATH $NB_INDIVIDUS $MIN $MAX $HEATMAP $FASTA #transfert des choix de l'utilisateur au script R
  
else
  echo "Hierarchique selectionne" ; echo ""
  
  echo "Vous avez le choix parmi ces methodes : "
    echo "-\"ward.D\" "
    echo "-\"ward.D2\""
    echo "-\"complete\""
    echo "-\"single\""
    echo "-\"average (UPGMA)\""
    echo "-\"mcquitty\""
    echo "-\"median\""
    echo "-\"centroid\""
  ###TMP### read -p "Quelle methode utiliser ? " METHODE ; echo ""
  ###TMP### read -p "Combien de groupes voulez-vous ? " NB_GROUPES ; echo ""
  
  METHODE="ward.D"
  NB_GROUPES=10
  
  if [ "y" == $FASTA ]
  then
    # read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " DATA_FASTA
    FASTA="403_seq.fasta"
  
    while [ "n" != $FASTA ] && [ ! -f $FASTA_PATH ]
    do
      echo "Fichier non trouve dans $DATA_DIR (taper 'n' pour ne pas avoir les fichiers fasta)"
      # read -p "Veuillez donner le nom du fichier fasta repertoriant vos individus : " DATA_FASTA
      FASTA="403_seq.fasta"
    done
  fi
  
  Rscript $MAIN_R $FIC_PATH $NB_INDIVIDUS $MIN $MAX $HEATMAP $FASTA $METHODE $NB_GROUPES #transfert des choix de l'utilisateur au script R
fi












