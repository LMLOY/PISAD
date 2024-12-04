#!/usr/bin/python
from __future__ import division
import os
import random
import math
from random import choice
import defSeq as dS
import re
import singleParameter as dI

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


owd= os.getcwd()

def aaList():
    os.chdir(owd)
    #os.chdir(home_dir)
    aa_list=[]
    with open("AllSeqData.txt", "r") as AllSequence:
        line=[lines for lines in AllSequence.readlines() if lines.strip()]
        for x in line:
            aa_list+=[x.split()[0]]
    return aa_list

def seq(seq_length, aaWV):
    sequence = "".join(random.choices(population=dS.ParAll, weights=aaWV, k=seq_length))
    return sequence

def PhilicRes(sequence):
    PhilicTotal=0
    for x in dS.ParPhilic:
        #PhilicTotal=PhilicTotal
        PhilicResid = sequence.count(x)
        PhilicTotal=PhilicResid+PhilicTotal
    print("Philic Residue in " + str(sequence) + ": " + str(PhilicTotal)) 
    return PhilicTotal

def PhobicRes(sequence):
    PhobicTotal=0
    for x in dS.ParPhobic:
        PhobicTotal=PhobicTotal
        PhobicResid = sequence.count(x)
        PhobicTotal=PhobicResid+PhobicTotal
    print("Phobic Residue in " + str(sequence) + ": " + str(PhobicTotal))
    return PhobicTotal

def changePhilicMin(seq_length, sequence, aaWV):
    minimum_philic = int(minper_philic*seq_length)
    changedPhilicMin = minimum_philic-PhilicRes(sequence)
    XPhilicList=[]
    for X in range(len(str(sequence))):
        if sequence[X] in dS.ParPhilic:
            pass
        else:
            XPhilicList+=[X]
    print("PhilicList: "+str(XPhilicList))
    sequence, indexAA = changingseq(sequence, seq_length, changedPhilicMin, XPhilicList, "IncPhilic", aaWV)
    print("after PhilicMin: "+str(sequence))
    return sequence
    
def changePhilicMax(seq_length, sequence, aaWV):
    if maxper_philic == []:
        sequence = sequence
    else: 
        maximum_philic = int(maxper_philic*seq_length)
        changedPhilicMax = PhilicRes(sequence)-maximum_philic
        XPhilicList=[]
        for X in range(len(str(sequence))):
            if sequence[X] in dS.ParPhilic:
                XPhilicList+=[X]
        print("PhilicList: "+str(XPhilicList))
        #indexPhilicList=[]
        sequence, indexAA = changingseq(sequence, seq_length, changedPhilicMax, XPhilicList, "RedPhilic", aaWV)
        PhilicResidue = PhilicRes(sequence)
        print("after PhilicMax: "+str(sequence))
    return sequence

def changePhobic(sequence, seq_length, aaWV):
    changedPhobicNum= PhobicRes(sequence) - int(maxper_phobic*seq_length)
    maximum_philic=int(maxper_philic*seq_length)
    changePhilicNum = maximum_philic-PhilicRes(sequence)
    sequence = sequence
    if maxper_philic == []:
        XPhobicList=[]
        for X in range(len(str(sequence))):
            if sequence[X] in dS.ParPhobic:
                XPhobicList+=[X]
        print("PhobicList: "+str(XPhobicList))
        sequence, indexAA = changingseq(sequence, seq_length, changedPhobicNum, XPhobicList, "RedPhobic", aaWV)
    else:
        PhilicResidue = PhilicRes(sequence)
        if PhilicResidue == maximum_philic:
            XPhobicList=[]
            for X in range(len(str(sequence))):
                if sequence[X] in dS.ParPhobic:
                    XPhobicList+=[X]
            sequence, indexAA = changingseq(sequence, seq_length, changedPhobicNum, XPhobicList, "Neutral", aaWV)        
        if PhilicResidue < maximum_philic:
            if changedPhobicNum > changePhilicNum:
                XPhobicList=[]
                for X in range(len(str(sequence))):
                    if sequence[X] in dS.ParPhobic and X != 0:
                        XPhobicList+=[X]
                sequence, indexAA = changingseq(sequence, seq_length, changePhilicNum, XPhobicList, "RedPhobic", aaWV)
                changedLeft=changedPhobicNum-changePhilicNum
                XPhobicList=[]
                for X in range(len(str(sequence))):
                    if sequence[X] in dS.ParPhobic:
                        XPhobicList+=[X]
                sequence, indexAA = changingseq(sequence, seq_length, changedLeft, XPhobicList, "Neutral", aaWV)
            else:
                XPhobicList=[]
                for X in range(len(str(sequence))):
                    if sequence[X] in dS.ParPhobic and X != 0:
                        XPhobicList+=[X]
                sequence, indexAA = changingseq(sequence, seq_length, changedPhobicNum, XPhobicList, "RedPhobic", aaWV)
    print("after Phobic: "+str(sequence))
    return sequence

def changingseq(sequence, seq_length, par1, par2, par3, aaWV):
    indexChangedList=[]
    for j in range (par1):
        indexChanged=choice(par2)
        while indexChanged in indexChangedList:
            indexChanged=choice(par2)
        indexChangedList+=[indexChanged]
    sortChangedList=sorted(indexChangedList)  
    print(sortChangedList)
            
    par3List = ["IncPhilic", "RedPhilic", "RedPhobic", "Neutral"]
    par4Helix = [dS.PhilicHelix, dS.NegPhilicHelix, dS.NegPhobicHelix, dS.NeutralHelix]
    par4Beta = [dS.PhilicBeta, dS.NegPhilicBeta, dS.NegPhobicBeta, dS.NeutralBeta]
    par4Turn = [dS.PhilicTurn, dS.NegPhilicTurn, dS.NegPhobicTurn, dS.NeutralTurn]
    par4Random = [dS.ParPhilic, dS.ParNoPhilic, dS.ParNoPhobic, dS.ParNeutral]
     
    if (seq_length-1) in sortChangedList: 
        for j in range (len(sortChangedList)-1):
            if par3 in par3List:
                par4 = []
                if sequence[sortChangedList[j]] in dS.ParHelix:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Helix[a]
                if sequence[sortChangedList[j]] in dS.ParBeta:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Beta[a]
                if sequence[sortChangedList[j]] in dS.ParTurn:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Turn[a]
                if par4 == []:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Random[a]
                else:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Random[a]
            if par3 not in par3List:
                par4 = par3
            aaWV1=[]
            for a in range(len(dS.ParAll)):
                if dS.ParAll[a] in par4:
                    aaWV1+=[aaWV[a]]
            sequence=sequence[:sortChangedList[j]]+"".join(random.choices(population=par4, weights=aaWV1, k=1))+sequence[sortChangedList[j]+1:]
        if par3 in par3List:
            par4 = []
            if sequence[sortChangedList[j]] in dS.ParHelix:
                for a in range(len(par3List)):
                    if par3 == par3List[a]:
                        par4 = par4Helix[a]
            if sequence[sortChangedList[j]] in dS.ParBeta:
                for a in range(len(par3List)):
                    if par3 == par3List[a]:
                        par4 = par4Beta[a]
            if sequence[sortChangedList[j]] in dS.ParTurn:
                for a in range(len(par3List)):
                    if par3 == par3List[a]:
                        par4 = par4Turn[a]
            if par4 == []:
                for a in range(len(par3List)):
                    if par3 == par3List[a]:
                        par4 = par4Random[a]
            else:
                for a in range(len(par3List)):
                    if par3 == par3List[a]:
                        par4 = par4Random[a]
        if par3 not in par3List:
            par4 = par3
        aaWV1=[]
        for k in range(len(par4)): 
            for a in range(len(dS.ParAll)):
                if dS.ParAll[a] == par4[k]:
                    aaWV1+=[aaWV[a]]
        sequence=sequence[:(seq_length-1)]+"".join(random.choices(population=par4, weights=aaWV1, k=1))

    else:
        for j in range (len(sortChangedList)):
            if par3 in par3List:
                par4 = []
                if sequence[sortChangedList[j]] in dS.ParHelix:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Helix[a]
                if sequence[sortChangedList[j]] in dS.ParBeta:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Beta[a]
                if sequence[sortChangedList[j]] in dS.ParTurn:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Turn[a]
                if par4 == []:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Random[a]
                else:
                    for a in range(len(par3List)):
                        if par3 == par3List[a]:
                            par4 = par4Random[a]
            if par3 not in par3List:
                par4 = par3
            aaWV1=[]
            for a in range(len(dS.ParAll)):
                if dS.ParAll[a] in par4:
                    aaWV1+=[aaWV[a]]
            sequence=sequence[:sortChangedList[j]]+"".join(random.choices(population=par4, weights=aaWV1, k=1))+sequence[sortChangedList[j]+1:]
    return sequence, sortChangedList

def NGPSRes(sequence):
    NGPSTotal=0
    for x in ["N","G","P","S"]:
        NGPSTotal=NGPSTotal
        NGPSResid = sequence.count(x)
        NGPSTotal=NGPSResid+NGPSTotal
    return NGPSTotal
    
def RKRes(sequence):
    RKTotal=0
    for x in ["R","K"]:
        RKTotal=RKTotal
        RKResid = sequence.count(x)
        RKTotal=RKResid+RKTotal
    return RKTotal
    
def DERes(sequence):
    DETotal=0
    for x in ["D","E"]:
        DETotal=DETotal
        DEResid = sequence.count(x)
        DETotal=DEResid+DETotal
    return DETotal

def CV(seq_length, parCV1, parCV2, parCV3):
    aaCV=((15.43*parCV1)/seq_length)-29.56*(abs(((parCV2-parCV3)/seq_length)-0.03))
    detCV=aaCV-1.71
    probCV=0.4934+(0.276*abs(detCV))-(0.0392*((detCV)**2))
    if detCV > 0:
        binCV=0
    else:
        binCV=1
    print("binary: "+str(binCV))
    return binCV, detCV

def II(seq_length, parII):
    weightTotal=0
    for l in range(len(str(parII))-1):
        weightRes=parII[l:l+2]
        weightTotal=weightTotal
        #print(weightTotal)
        #print(weightRes)
        os.chdir(owd)
        #os.chdir(home_dir)
        with open("InstabilityWeight.txt", "r") as weightFile:
            for index, line in enumerate(weightFile, 1):
                if weightRes in line:
                    indexWeight=index
                    #print(indexWeight)
                    weightValue=float(line.split()[1])
                    #print("weightvalue:" +str(weightValue))
        weightTotal=weightValue+weightTotal
    print("Total weight value:" +str(weightTotal))
    InstabilityIndex=10*weightTotal/seq_length
    if InstabilityIndex > 40:
        binII=0
    else:         
        binII=1
    print("binary: "+str(binII))
    return binII, InstabilityIndex

def duplicate(sequence, letter):
    residue = re.finditer(letter, sequence)
    position = [resid.start() for resid in residue]
    return position






    
