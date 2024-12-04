#!/usr/bin/python
from __future__ import division
import os
import random
from random import choice
import math
import defSeq as dS
import defsingleEquation as dE
import defsingleLogic as dL
import re
import singleParameter as dI
import defsingleMut as dM
import sys

#Getting Parameter from singleParameter.py
home_dir = dI.home_dir
OutSeq_dir = dI.OutSeq_dir
psipred_dir = dI.psipred_dir
Mut_dir = dI.Mut_dir
seq_length = dI.seq_length
HHCriteria = dI.HHCriteria
minper_philic=dI.minper_philic
maxper_philic=dI.maxper_philic
maxper_phobic=dI.maxper_phobic
SSCriteria = dI.SSCriteria
Perc_Helix = dI.Perc_Helix
Perc_Beta = dI.Perc_Beta
Perc_Turn = dI.Perc_Turn
aaWV = dI.aaWV

mut_file = sys.argv[1]
Consecutive = dI.Consecutive
mutPosition = dI.mutPosition
mutWay = dI.mutWay
mutRes = dI.mutRes
    
owd= os.getcwd()
os.chdir(owd)

def SeqCall(mut_file):
    os.chdir(Mut_dir)
    seqFile=open(mut_file)
    seqLine=seqFile.readlines()
    seq=seqLine[1]
    os.chdir(owd)
    return seq

aa_list=dE.aaList()

#Calling Sequence
aaSequence = SeqCall(mut_file)
oriSequence = SeqCall(mut_file)

if SSCriteria == "Yes" or SSCriteria == "yes":
    for x in range(len(dS.ParAll)):
        if dS.ParAll[x] in dS.ParHelix:
            aaWV[x] = aaWV[x]+Perc_Helix
        if dS.ParAll[x] in dS.ParBeta:
            aaWV[x] = aaWV[x]+Perc_Beta
        if dS.ParAll[x] in dS.ParTurn:
            aaWV[x] = aaWV[x]+Perc_Turn
        else:
            aaWV[x] = aaWV[x]
if SSCriteria == "No" or SSCriteria == "no":
    aaWV = aaWV

#Normalization of Weight Value
WVTotal=0
for x in range(len(aaWV)):
    WVTotal = aaWV[x]+WVTotal
for x in range(len(aaWV)):
    aaWV[x] = aaWV[x]/WVTotal

if Consecutive not in ["Yes", "yes", "No", "no"]:
    raise ValueError("Error in Consecutive Parameter. True or False only.")
if len(mutPosition) > mutRes:
    raise ValueError("The selected mutated position is larger than the set number for mutation")
if mutRes > seq_length-1 or mutRes < 0:
    raise ValueError("number of mutated letters larger than sequence seq_length or is a negative value")            
if mutRes == 0:
    raise ValueError("The number of residue to be mutated is 0. Please insert new mutRes value in InputParameter.py.")


if mutWay == 0:
    mutWay=choice([1, 2])

aaSequence = dM.mutCases(aaSequence, seq_length, mutWay, mutRes, mutPosition, aaWV)
aaSequence, aa_list, binCV, binII = dM.mutCasesEval(aaSequence, oriSequence, aa_list, seq_length, aaWV)

while binCV==0 or binII==0:
    aaSequence=oriSequence
    mutWay=dI.mutWay
    if mutWay == 0:
        mutWay=choice([1, 2])
    aaSequence = dM.mutCases(aaSequence, seq_length, mutWay, mutRes, mutPosition, aaWV)
    aaSequence, aa_list, binCV, binII = dM.mutCasesEval(aaSequence, oriSequence, aa_list, seq_length, aaWV)
    print("after failed:"+str(aaSequence))

if binCV==1 and binII==1:
    os.chdir(psipred_dir)
    f= open("aa_temp.fasta","w+")
    f.write(">aaTemp \n")
    f.write("".join(aaSequence))
    os.chdir(owd)
    os.chdir(home_dir)
    with open("AllSeqData.txt", "a") as file:
        file.write(str(aaSequence)+"\t")