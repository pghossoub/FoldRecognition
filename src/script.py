#!/usr/bin/env python3
# -*- coding:Latin-1 -*-

import sys
import argparse
from shutil import copy
import subprocess
import os

def get_parametres():

    parametres = {"OrionPath": "/home/sdv/m2bi/pghossoub/Download/", "FirstLine": 9, "LastLine" : 20}

    parser = argparse.ArgumentParser()

    parser.add_argument("--op", help="Chemin du répertoire contenant le répertoire orion_files.")
    parser.add_argument("-f", help="Ligne de début de lecture du fichier posAlignUp25-35.txt contenant les groupes de familles d'alignement positifs")
    parser.add_argument("-l", help="Ligne de fin de lecture du fichier posAlignUp25-35.txt contenant les groupes de familles d'alignement positifs")



    args = parser.parse_args()

    return parametres

def read_align(first, last):

	"""
	Lit les noms des alignements positifs contenus dans posAlignUp25-35.txt 
	qui sont les noms des fichiers souhaités pour notre learning de la base
	de donnée orion.
	Stocke ces noms dans une liste de liste:[groupe][famille] 
	"""

	filePATH = "data/posAlignUp25-35.txt"
	#listGroup = []
	listFamille = []

	with open(filePATH, 'r') as fin:
		#for line in fin:
		lines = fin.readlines()
		for i in xrange(int(first),int(last)):
			#listFamille.append(line.split('\t')[:-1])
			listFamille.append(lines[i].split('\t')[:-1])


	#print listGroup
	return listFamille



def align(orionPath, listFamille):
	"""
	Cherche dans la base de donnée orion (en local) les noms des
	alignements.fas de la listFamille considérée.

	Les aligne avec mafft avec une nouvelle séquence d'un autre
	groupe pour générer un alignement négatif.

	Supprime le début du fichier out de l'alignement mafft
	pour enlever la partie positive.
	"""

	command = "src/mafft-7.313-without-extensions/binaries/ginsi"


	for i in xrange(len(listFamille)):
		for j in xrange(len(listFamille[i])):

			pos = listFamille[i][j]
			if(i < len(listFamille)-1):
				neg = listFamille[i+1][j]
			else:
				neg = listFamille[i-1][j]

			pathAlignPos= orionPath + "orion_files/upload_orion_files/" + pos + "/" + pos + "_mafft2.fas"
			pathAlignNeg= orionPath + "orion_files/upload_orion_files/" + neg + "/" + neg + "_mafft2.fas"



			

			#CMD_list = ["ls", "-l"]#arg 
			#subprocess.call(CMD_list)

			CMD_list = [command, "--thread", "-1", "--addfull", pathAlignNeg, "--keeplength", pathAlignPos, ">", "data/bd/align_negatif/"+ pos + "NEG.mafft2.fas"]
			if(os.path.isfile(pathAlignNeg) and os.path.isfile(pathAlignPos)):

				subprocess.call(CMD_list)

				copy(pathAlignPos, "data/bd/align_positif")



if __name__ == "__main__":

	parametres = get_parametres()
	#print parametres
	#print parametres["FirstLine"]

	listFamille = read_align(parametres["FirstLine"], parametres["LastLine"])
	print listFamille

	#print parametres["OrionPath"]
	align(parametres["OrionPath"], listFamille)
