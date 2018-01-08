#!/usr/bin/env python3
# -*- coding:Latin-1 -*-


#def extractAAINDEX():
if __name__ == "__main__":

    listRepresentAAINDEX = []

    with open("data/AAINDEX58", "r") as fin:
        #lines = fin.readlines()
        for line in fin:
            listRepresentAAINDEX.append(line.split("\n")[:-1])


    listMots = []

    with open("data/selected_aaindex1_reformated", "r") as fin:
        #lines = fin.readlines()
        for line in fin:
            listMots.append(line.split())

    toWrite =""
    with open("selected_AAINDEX58_reformated", "w") as fout:
        for i in xrange(len(listMots)):
            if(listMots[i][0] in str(listRepresentAAINDEX) or i == 0):
                print "trouvé"
                toWrite = " ".join(listMots[i])
                fout.write(toWrite + "\n")
                toWrite =""


#if __name__ == "__main__":

#    extractAAINDEX()