#!/bin/bash

ROOT=$(pwd)
DATA_DIR="$ROOT/data/"
SCRIPT="/scripts/"
SCRIPT_DIR="$ROOT$SCRIPT"
VLD="/decode_src/"
VLD_DIR="$ROOT$VLD"

read -p "Avez-vous deja une matrice de distances (fichier .raw) ? (y/n)" SHORTCUT
if [ "n" == $SHORTCUT ]
then
  #Fichier fasta en entree
  read -p "Donnez le nom de votre fichier fasta : " FASTA_FIC_NAME
  FASTA_FIC_PATH="$DATA_DIR$FASTA_FIC_NAME"
  
  #Application de VLD sur ce fichier
  FIC=$(basename $FASTA_FIC_PATH)
  FIC_NAME="${FIC%.*}_VLD.fasta"
  
  #Creation fichier fasta VLD
  echo "Reecriture de l'alphabet de $FASTA_FIC_NAME avec VLD"
  VLD_FASTA_FIC_PATH="$DATA_DIR$FIC_NAME"
  ."$VLD"decode $FASTA_FIC_PATH $VLD_FASTA_FIC_PATH
  
  #Creation matrice de distances
  echo "Creation de la matrice de distances"
  VLD_FIC=$(basename $VLD_FASTA_FIC_PATH)
  RAW_VLD_FIC_NAME="${VLD_FIC%.*}.raw"
  RAW_VLD_FIC_PATH="$DATA_DIR$RAW_VLD_FIC_NAME"
  ."$VLD"distauto -f r -s p $VLD_FASTA_FIC_PATH $RAW_VLD_FIC_PATH
  
  echo "Placez maintenant la derniere ligne du fichier $RAW_VLD_FIC_PATH en premier"
  read -p "Entree pour continuer" REP
fi

echo "" ; echo ""
#Application de la classification ensuite
."$SCRIPT"clustering.sh
