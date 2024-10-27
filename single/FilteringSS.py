#!/usr/bin/python
from __future__ import division
import os
import sys
import InputParameter as dI 
import math
import re
import defEquation as dE

#Filtering parameter
#1. total number of coil in structure is no more than 30% of total sequence length or total chain length (for polypeptide)
#2. number of coil in a row (e.g. CCCC) no more than 30% of total sequence length or total chain length (for polypeptide)

seq_dir=dI.seq_dir
home_dir=dI.home_dir
psipred_dir=dI.psipred_dir
polypeptides = dI.polypeptides
seq_length = dI.seq_length

num_chain = dI.num_chain
length_seq1 = dI.length_seq1
length_linker1 = dI.length_linker1
length_seq2 = dI.length_seq2
length_linker2 = dI.length_linker2
length_seq3 = dI.length_seq3
type_linker1 = dI.type_linker1
type_linker2 = dI.type_linker2

coilPercentage = dI.coilPercentage
coilFlex = dI.coilFlex
coilRigid = dI.coilRigid
i=sys.argv[5]
typeIMC=str(sys.argv[6])
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
if polypeptides == "False":
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
            #print(Coilposition[x])
            if Coilposition[x+1] == Coilposition[x]+1:
                coilCount+=1
            else:
                coilCount=1
            coilCountList+=[coilCount]
            #print(coilCount)
            #print(coilCountList)
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
                j = sys.argv[7]
                saveas = str(i)+"_mutation"+str(j)
            if typeIMC == "Crossover" :
                saveas = "crossover"+str(i)
            with open("SequenceFile.txt", "a") as file:
                file.write(str(saveas)+"\t"+str(aa)+"\t"+str(charge)+"\n")
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

length_list=[]
aaChain=[]
ssChain=[]
charge1=""
charge2=""
                         
if polypeptides == "True":
    length_list=[length_seq1, length_linker1, length_seq2]
    if num_chain == 3:
        length_list+=[length_linker2, length_seq3]
    for x in range(len(length_list)):
        if x in [0, 2, 4]:
            maxCoil = int(coilPercentage*length_list[x])
            maxCoilChain+=[maxcoil]
            if x == 0:
                aaChain+=[aa[:length_seq1]]
                ssChain+=[ss[:length_seq1]]
            if x == 2:
                aaChain+=[aa[length_seq1+length_linker1:length_seq1+length_linker1+length_seq2]]
                ssChain+=[ss[length_seq1+length_linker1:length_seq1+length_linker1+length_seq2]]
            if x == 4:
                aaChain+=[aa[length_seq1+length_linker1+length_seq2+length_linker2:]]
                ssChain+=[ss[length_seq1+length_linker1+length_seq2+length_linker2:]]
        if x == 1:
            if type_linker1 == "flexible":
                maxCoil = int(coilFlex * length_list[x])
            if type_linker1 == "rigid":
                maxCoil = int(coilRigid * length_list[x])
            if type_linker1 == "random":
                maxCoil = length_linker1
            maxCoilChain+=[maxCoil]
            aaChain+=[aa[length_seq1:length_seq1+length_linker1]]
            ssChain+=[ss[length_seq1:length_seq1+length_linker1]]
        if x == 3:
            if type_linker2 == "flexible":
                maxCoil = int(coilFlex * length_list[x])
            if type_linker2 == "rigid":
                maxCoil = int(coilRigid * length_list[x])
            if type_linker2 == "random":
                maxCoil = length_linker2
            maxCoilChain+=[maxCoil]
            aaChain+=[aa[length_seq1+length_linker1+length_seq2:length_seq1+length_linker1+length_seq2+length_linker2]]
            ssChain+=[ss[length_seq1+length_linker1+length_seq2:length_seq1+length_linker1+length_seq2+length_linker2]]             
    ResultList=[]
    for x in range(len(aaChain)):
        Coil = re.finditer('C', ssChain[x])
        Coilposition = [Coils.start() for Coils in Coil]
        if maxCoilChain[x] > 10 :
            if length_list[x] <30:
                 maxCoilChain[x]=5
            if length_list[x] > 29 and length_list[x] < 71:
                maxCoilchain[x]=9
            if length_list[x]>70:
                maxCoilChain[x]=21

        if len(Coilposition) > maxCoilChain[x]:
            ResultList+=["F"]
        else:
            coilCount=1
            coilCountList=[]
            for y in range(len(Coilposition)-1):
                #print(Coilposition[x])
                if Coilposition[y+1] == Coilposition[y]+1:
                    coilCount+=1
                else:
                    coilCount=1
                coilCountList+=[coilCount]
                #print(coilCount)
                #print(coilCountList)
            if maxCoilChain[x]+1 in coilCountList:
                ResultList+=["F"]
            else:
                ResultList+=["P"]
    if "F" in ResultList:
        Result="F"
    if "F" not in ResultList:
        Result="P"
    os.chdir(owd)
    os.chdir(home_dir)
    with open("AllSeqData.txt", "r") as file:
        line=[lines for lines in file.readlines() if lines.strip()]
        for x in range(len(line)):
            #print(str(line[x].split()[0]))
            if str(line[x].split()[0]) == str(aa):
                charge=line[x].split()[1]
                charge1=line[x].split()[2]
                if num_chain == 3:
                    charge2=line[x].split()[3]
                line[x]=line[x].replace('\n', '\t')+(str(Result)+ "\n")
    with open("AllSeqData.txt", "w") as file:
        file.writelines(line)
    if Result == "P":
        if typeIMC == "Initial":
            saveas = "aa"+str(i)
        if typeIMC == "Mutation":
            j = sys.argv[7]
            saveas = str(i)+"_mutation"+str(j)
        if typeIMC == "Crossover" :
            saveas = "crossover"+str(i)
        with open("SequenceFile.txt", "a") as file:
            if num_chain == 2:
                file.write(str(saveas)+"\t"+str(aa)+"\t"+str(charge)+"\t"+str(charge1)+"\n")
            if num_chain == 3:
                file.write(str(saveas)+"\t"+str(aa)+"\t"+str(charge)+"\t"+str(charge1)+"\t"+str(charge2)+"\n")
        os.chdir(owd)
        os.chdir(seq_dir)
        f = open(str(saveas)+".fasta", "w+")
        f.write(str(saveas)+"\n")
        f.write(str(aa))
        f.close()
