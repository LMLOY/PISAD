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
import defsingleCross as dC

owd= os.getcwd()

#Getting Parameter from singleParameter.py
home_dir = dI.home_dir
OutSeq_dir = dI.OutSeq_dir
psipred_dir = dI.psipred_dir
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
Cross_dir = dI.Cross_dir

aa_list=dE.aaList()

def SeqCall(name):
    os.chdir(Cross_dir)
    seqFile=open(str(name)+".fasta")
    seqLine=seqFile.readlines()
    seq=seqLine[1]
    os.chdir(owd)
    return seq

os.chdir(Cross_dir)
crossfile=[]
for files in os.listdir("."):
    filename=files.rsplit(".")[0]
    crossfile+=[filename]
print(crossfile)
os.chdir(owd)
file1=choice(crossfile)
file2=choice(crossfile)
while file2 == file1:
    file2=choice(crossfile)
sequence1=SeqCall(file1)
sequence2=SeqCall(file2)
print(sequence1)
print(sequence2)
sequence1, sequence2, sameAAList, redo = dC.evaluateFindSeq(sequence1, sequence2, seq_length)
while redo == True:
    file1=choice(crossfile)
    file2=choice(crossfile)
    while file2 == file1:
        file2=choice(crossfile)
    sequence1=SeqCall(file1)
    sequence2=SeqCall(file2)
    sequence1, sequence2, sameAAList, redo = dC.evaluateFindSeq(sequence1, sequence2, seq_length)
crossSequence, aa_list = dC.crossProcess(sequence1, sequence2, sameAAList, aa_list, seq_length)
PhilicResidue, PhobicResidue, binCV, CanonicalValue, binII, InstabilityIn = dC.crossEval(crossSequence, seq_length)

while binCV == 0 or binII == 0:
    print("restart everything")
    file1=choice(crossfile)
    file2=choice(crossfile)
    while file2 == file1:
        file2=choice(crossfile)
    sequence1=SeqCall(file1)
    sequence2=SeqCall(file2)
    sequence1, sequence2, sameAAList, redo = dC.evaluateFindSeq(sequence1, sequence2, seq_length)
    while redo == True:
        file1=choice(crossfile)
        file2=choice(crossfile)
        while file2 == file1:
            file2=choice(crossfile)
        sequence1=SeqCall(file1)
        sequence2=SeqCall(file2)
        sequence1, sequence2, sameAAList, redo = dC.evaluateFindSeq(sequence1, sequence2, seq_length)
    crossSequence, aa_list = dC.crossProcess(sequence1, sequence2, sameAAList, aa_list, seq_length)
    PhilicResidue, PhobicResidue, binCV, CanonicalValue, binII, InstabilityIn = dC.crossEval(crossSequence, seq_length)
    
if binCV == 1 and binII == 1:
    os.chdir(psipred_dir)
    f= open("aa_temp.fasta","w+")
    f.write(">aaTemp \n")
    f.write("".join(crossSequence))
    os.chdir(owd)
    os.chdir(home_dir)
    with open("AllSeqData.txt", "a") as file:
        file.write(str(crossSequence)+"\t")
