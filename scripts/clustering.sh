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

read -p "Veuillez indiquer le nom de votre matrice de distances (fichier .raw) : " FIC_NAME
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
  echo "Veuillez placer votre fichier dans le dossier data avant de relancer le programme"
  exit
else
  echo "Fichier de donnees existant dans $DATA_DIR !"
  echo ""
fi

#A ce stade, le fichier a ete trouve





















