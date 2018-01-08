#!/usr/bin/env python3
# -*- coding:Latin-1 -*-

import sys
import argparse
import shutil
import subprocess
import os

def get_parametres():

    parametres = {"OrionPath": "/home/sdv/m2bi/pghossoub/Download/", "FirstLine": 9, "LastLine" : 10}

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
        for i in xrange(int(first),int(last)+1): #+1 car la dernière ligne donnée serait exclue
            #listFamille.append(line.split('\t')[:-1])
            listFamille.append(lines[i].split('\t')[:-1])


    #print listGroup
    return listFamille



def false_align(orionPath, listFamille):
    """
    Cherche dans la base de donnée orion (en local, au chemin OrionPath) 
    les noms des alignements.fas de la listFamille considérée.

    Les aligne avec mafft avec une nouvelle séquence d'un autre
    groupe pour générer un alignement négatif.

    Supprime le début du fichier out de l'alignement mafft
    pour enlever la partie positive.
    """

    print "\nfalse alignement\n"
    command = "src/mafft-7.313-without-extensions/binaries/ginsi"


    for i in xrange(len(listFamille)):
        #for j in xrange(len(listFamille[i])):
        for j in xrange(len(listFamille[i])):

            pos = listFamille[i][j]
            print pos

            if(i < len(listFamille)-1):
                neg = listFamille[i+1][0] #j = 0 par exemple, c'est un d'un autre groupe
            else:
                neg = listFamille[i-1][0] #j = 0 par exemple, c'est un d'un autre groupe

            print neg

            pathAlignPos= orionPath + "orion_files/upload_orion_files/" + pos + "/" + pos + "_mafft2.fas"
            pathAlignNeg= orionPath + "orion_files/upload_orion_files/" + neg + "/" + neg + "_mafft2.fas"
            print pathAlignPos
            print pathAlignNeg

            #CMD_list = ["ls", "-l"]#arg 
            #subprocess.call(CMD_list)

            #CMD_list = [command, "--thread", "-1", "--addfull", pathAlignNeg, "--keeplength", pathAlignPos, "> data/bd/align_negatif/"+ pos + "NEG.mafft2.fas"]
            #CMD_list = [command, "--thread", "-1", "--addfull", pathAlignNeg, "--keeplength", pathAlignPos] #marche, pas de out

            sortieTMP = "data/bd/" + pos + "TMP.mafft2.fas"

            CMD_list = command + " --thread -1 --addfull "+ pathAlignNeg + " --keeplength " + pathAlignPos + " > " + sortieTMP




            if(os.path.isfile(pathAlignNeg) and os.path.isfile(pathAlignPos)):
            	print "\nmafft\n"
                #sortie = "data/bd/align_negatif/"+ pos + "NEG.mafft2.fas"
                #subprocess.call(CMD_list, stdout =  sortie)

                subprocess.call(CMD_list, shell = True)
                shutil.copy(pathAlignPos, "data/bd/align_positif")

                keepOnlyNegative("data/bd/align_positif/" + pos + "_mafft2.fas" , sortieTMP, pos)




def compte_ligne(fin):

    n = 0
    for line in fin:
        n += 1

    return n


def keepOnlyNegative(nomFichierPos, nomFichierTmp, nomAlignement):

    newFichier = "data/bd/align_negatif/" + nomAlignement + "NEG_mafft2.fas"

    with open(nomFichierPos, 'r') as fin:
        toRemove = compte_ligne(fin)


    with open(nomFichierTmp, 'r') as fin:
        with open(newFichier, 'w') as fout:
            for i, line in enumerate(fin):
                if i >= toRemove:
                    fout.write(line)
    


def align_To_Profile_AAINDEX():

    pathPos = "data/bd/align_positif/"
    pathNeg = "data/bd/align_negatif/"

    fasta2vector = "perl src/gelly/fasta2vector_wgap.pl"
    aaindex58 = "data/selected_AAINDEX58_reformated"

    fileInDir = os.listdir(pathPos)
    for f in fileInDir:
        if os.path.isfile(os.path.join(pathPos, f)):
            out =  pathPos + f + ".AAINDEX58"
            CMD_list = fasta2vector + " " + pathPos + f + " " + aaindex58 + " > " + out
            print CMD_list
            subprocess.call(CMD_list , shell = True)

    fileInDir = os.listdir(pathNeg)
    for f in fileInDir:
        if os.path.isfile(os.path.join(pathNeg, f)):
            out =  pathNeg + f + ".AAINDEX58"
            CMD_list = fasta2vector + " " + pathNeg + f + " " + aaindex58 + " > " + out
            print CMD_list
            subprocess.call(CMD_list , shell = True)



if __name__ == "__main__":

    
    parametres = get_parametres()
    #print parametres
    #print parametres["FirstLine"]

    listFamille = read_align(parametres["FirstLine"], parametres["LastLine"])
    print listFamille

    #print parametres["OrionPath"]
    false_align(parametres["OrionPath"], listFamille)
    
    align_To_Profile_AAINDEX()