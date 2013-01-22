#!/usr/bin/env python
"""
Created on Mon Oct  8 10:13:40 2012

@author: jim
"""

import sys
import argparse
import subprocess
import csv
import re
import random
import os

import Bio
from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers, and produces a relative abundance table.')


parser.add_argument('--nucs', type=str, dest='sNucs', help='Enter the path and name of the nucleotide file.')
parser.add_argument('--genomes', type=str, dest='sGenomes', help='Enter the path and name of the file containing the names of the full genomes.')
parser.add_argument('--var', type=int, dest='iUB', help='Enter the upper bound for the U(1,UB) draws for each gene.',default=1)
parser.add_argument('-N', type=int, dest='iN', help='Enter the number of genes to include in the file.',default =10)
parser.add_argument('--gold', type=str, dest='sGold', help='Enter the path and name of the gold standard file.',default="goldstandard.txt")



#parser.add_argument('--pct', type=float, dest='dPct', help='Enter the starting percentage, <= .30.', default=.20)
#parser.add_argument('--log', type=str, dest='sLog',help='Enter the name of the log file')
#parser.add_argument('--map', type=str, dest='sMap',help='Enter the name of the map file from usearch')


#parser.add_argument('--props', type=str, dest='sProps', help='Enter the path and name of the file with gene names and proportions.')

args = parser.parse_args()

dictNucs = {}
aaGoldStandard = []

fileGS = open(args.sGold,'w')

#Load in data from nuc file
for seq in SeqIO.parse(args.sNucs, "fasta"):
    seq.id = "USR_" + seq.id + "_END"
    dictNucs[seq.id] = seq

#Choose N nucs randomly
setPickedGenes = random.sample(set(dictNucs.keys()), args.iN)


astrPickedGenomes = [line.strip() for line in open(args.sGenomes)]

for strGeneName in astrPickedGenomes:
    aGeneData = [strGeneName, random.lognormvariate(6.25,1), "NA","NA"]
    aaGoldStandard.append(aGeneData)



#Create an array of the percentages for the simulated genes
for strGeneName in setPickedGenes:
    aGeneData = [strGeneName, random.lognormvariate(1,args.iUB), len(dictNucs[strGeneName])]
    aGeneData.append(aGeneData[1]*aGeneData[2])
    aaGoldStandard.append(aGeneData)

fileGS.write("GeneName"+ "\t" + "Count" + "\t" + "Length" + "\n")
for aLine in aaGoldStandard:
    fileGS.write(aLine[0]+ "\t" + str(aLine[1]) + "\t" +str(aLine[2]) + "\t" + str(aLine[3]) )
    fileGS.write("\n")

"""
astrGeneRepeats = []

for aLine in aaGoldStandard:
    astrCopies =  aLine[1]*[aLine[0]]
    for x in astrCopies:
        astrGeneRepeats.append(x)

#print astrGeneRepeats
random.shuffle(astrGeneRepeats)
#print astrGeneRepeats


for strGene in astrGeneRepeats:
    SeqIO.write(dictNucs[strGene], sys.stdout, "fasta")

iTotBases = sum(zip(*aaGoldStandard)[3])
iTotGenes = sum(zip(*aaGoldStandard)[2])

fileGS.write("The shuffled and expanded genome file for " + args.sNucs + " contains" + "\n")

fileGS.write("Input genes: " + str(len(setPickedGenes)) + "\n")
fileGS.write("Gene copies: " + str(iTotGenes) + "\n")
fileGS.write("Total bases: " + str(iTotBases) + "\n")
"""
fileGS.close()




