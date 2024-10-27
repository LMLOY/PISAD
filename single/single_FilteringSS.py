#!/usr/bin/python
from __future__ import division
import os
import sys
import singleParameter as dI
import math
import re

#Filtering parameter
#1. total number of coil in structure is no more than 30% of total sequence length or total chain length (for polypeptide)
#2. number of coil in a row (e.g. CCCC) no more than 30% of total sequence length or total chain length (for polypeptide)

#Getting Parameter from singleParameter.py
home_dir = dI.home_dir
OutSeq_dir = dI.OutSeq_dir
psipred_dir = dI.psipred_dir
seq_length = dI.seq_length
coilPercentage = dI.coilPercentage

typeIMC=str(sys.argv[1])
i=sys.argv[2]
j=""

owd= os.getcwd()
os.chdir(psipred_dir)
ss=""
aa=""
Result=""
saveas=""

with open("aa_temp.ss2", "r") as ssFile:
    line=[lines for lines in ssFile.readlines() if lines.strip()]
    #print(len(line))
    for x in range(len(line)):
        if x < 1:
            pass
        else:
            ss=str(ss)+line[x].split()[2]
            aa=str(aa)+line[x].split()[1]
#print(ss)
#print(aa)
print(i)

maxCoil = int(coilPercentage*seq_length)
Coil = re.finditer('C', ss)
Coilposition = [Coils.start() for Coils in Coil]
if len(Coilposition) > maxCoil:
    Result="F"
    os.chdir(owd)
    os.chdir(home_dir)
    with open("AllSeqData.txt", "r") as file:
        line=[lines for lines in file.readlines() if lines.strip()]
        for x in range(len(line)):
            #print(str(line[x].split()[0]))
            if str(line[x].split()[0]) == str(aa):
                line[x]=line[x].replace('\n', '\t')+(str(ss)+"\t F \n")
    with open("AllSeqData.txt", "w") as file:
        file.writelines(line)
else:
    coilCount=1
    coilCountList=[]
    for x in range(len(Coilposition)-1):
        if Coilposition[x+1] == Coilposition[x]+1:
            coilCount+=1
        else:
            coilCount=1
        coilCountList+=[coilCount]

    if maxcoil > 10:
        if seq_length < 30:
            maxcoil1=5
        if seq_length > 29 and seq_length < 71:
            maxcoil1=9
        if seq_length > 70:
            maxcoil1=21
    if maxcoil1+1 in coilCountList:
        Result="F"
        os.chdir(owd)
        os.chdir(home_dir)
        with open("AllSeqData.txt", "r") as file:
            line=[lines for lines in file.readlines() if lines.strip()]
            for x in range(len(line)):
                print(str(line[x].split()[0]))
                if str(line[x].split()[0]) == str(aa):
                    line[x]=line[x].replace('\n', '\t')+(str(ss)+"\t F \n")
        with open("AllSeqData.txt", "w") as file:
            file.writelines(line)
    else:
        Result="P"
        os.chdir(owd)
        os.chdir(home_dir)
        with open("AllSeqData.txt", "r") as file:
            line=[lines for lines in file.readlines() if lines.strip()]
            for x in range(len(line)):
                print(str(line[x].split()[0]))
                if str(line[x].split()[0]) == str(aa):
                    charge=line[x].split()[1]
                    line[x]=line[x].replace('\n', '\t')+(str(ss)+"\t P \n")
        with open("AllSeqData.txt", "w") as file:
            file.writelines(line)
        if typeIMC == "Initial":
            saveas = "aa"+str(i)
        if typeIMC == "Mutation":
            j = sys.argv[3]
            saveas = str(i)+"_mutation"+str(j)
        if typeIMC == "Crossover" :
            saveas = "crossover"+str(i)
        with open("SequenceFile.txt", "a") as file:
            file.write(str(saveas)+"\t"+str(aa)+"\n")
        os.chdir(owd)
        os.chdir(seq_dir)
        f = open(str(saveas)+".fasta", "w+")
        f.write(str(saveas)+"\n")
        f.write(str(aa))
        f.close()
os.chdir(owd)
os.chdir(psipred_dir)
if os.path.exists("aa_temp.fasta"):
    os.remove("aa_temp.fasta")
if os.path.exists("aa_temp.ss2"):
    os.remove("aa_temp.ss2")
else:
    pass
