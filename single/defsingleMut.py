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


mutPosition = dI.mutPosition
mutWay = dI.mutWay
mutRes = dI.mutRes

owd= os.getcwd()
#os.chdir(seq_dir)

def mutCases(aaSequence, seq_length, mutWay, mutRes, mutPosition, aaWV):
    aaSequence = aaSequence
    print(aaSequence)
    if dI.Consecutive == "Yes" or dI.Consecutive == "yes":
        if mutPosition != []:
            mutatedList=[]
            for k in range(mutRes):
                mutatedList+=[mutPosition-1+k]
        else:
            XChanged_List=[]
            for x in range(len(str(aaSequence))-mutRes):
                XChanged_List+=[x]
            mutPosition=choice(XChanged_List)
            mutatedList=[]
            for k in range(mutRes):
                mutatedList+=[mutPosition+k]
    if dI.Consecutive == "No" or dI.Consecutive == "no":
        if mutPosition != []:
            mutatedList=mutPosition-1
        else:
            XChanged_List=[]
            for x in range(len(str(aaSequence))):
                XChanged_List+=[x]
            mutatedList=[]
            for k in range(mutRes):
                z=choice(XChanged_List)
                while z in mutatedList:
                    z=choice(XChanged_List)
                mutatedList+=[z]
    sortChangedList=sorted(mutatedList)
    print(sortChangedList)
    
    #Substitute with the same characteristic residue. 
    #Ex: if the selected residue in Neutral list then the substitute residue also from Neutral List
    if mutWay == 1:
        choiceResList=[]
        for k in range(len(sortChangedList)):
            print(k)
            if aaSequence[sortChangedList[k]] in dS.ParPhilic:
                if aaSequence[sortChangedList[k]] in dS.ParHelix:
                    choicePar = dS.PhilicHelix
                if aaSequence[sortChangedList[k]] in dS.ParBeta:
                    choicePar = dS.PhilicBeta
                if aaSequence[sortChangedList[k]] in dS.ParTurn:
                    choicePar = dS.PhilicTurn
                if aaSequence[sortChangedList[k]] not in [dS.ParHelix, dS.ParBeta, dS.ParTurn] or len(choicePar) < 2:
                    choicePar = dS.ParPhilic
            if aaSequence[sortChangedList[k]] in dS.ParPhobic:
                if aaSequence[sortChangedList[k]] in dS.ParHelix:
                    choicePar = dS.PhobicHelix
                if aaSequence[sortChangedList[k]] in dS.ParBeta:
                    choicePar = dS.PhobicBeta
                if aaSequence[sortChangedList[k]] in dS.ParTurn:
                    choicePar = dS.PhobicTurn
                if aaSequence[sortChangedList[k]] not in [dS.ParHelix, dS.ParBeta, dS.ParTurn] or len(choicePar) < 2:
                    choicePar = dS.ParPhobic
            if aaSequence[sortChangedList[k]] in dS.ParNeutral:
                if aaSequence[sortChangedList[k]] in dS.ParHelix:
                    choicePar = dS.NeutralHelix
                if aaSequence[sortChangedList[k]] in dS.ParBeta:
                    choicePar = dS.NeutralBeta
                if aaSequence[sortChangedList[k]] in dS.ParTurn:
                    choicePar = dS.NeutralTurn
                if aaSequence[sortChangedList[k]] not in [dS.ParHelix, dS.ParBeta, dS.ParTurn] or len(choicePar) < 2:
                    choicePar = dS.ParNeutral
                    
            aaWV1=[]
            for j in range(len(choicePar)): 
                for a in range(len(dS.ParAll)):
                    if dS.ParAll[a] == choicePar[j]:
                        aaWV1+=[aaWV[a]]
            print(choicePar)
            print(aaWV1)
            choiceRes = "".join(random.choices(population=choicePar, weights=aaWV1))
            while choiceRes == aaSequence[sortChangedList[k]]:
                choiceRes = "".join(random.choices(population=choicePar, weights=aaWV1))
            choiceResList+=[choiceRes]
            print(choiceRes)
        print(choiceResList)
        
        if seq_length-1 in sortChangedList:
            for j in range(len(sortChangedList)-1):
                aaSequence=aaSequence[:sortChangedList[j]]+str(choiceResList[j])+aaSequence[sortChangedList[j]+1:]
            aaSequence=aaSequence[:(seq_length-1)]+str(choiceResList[len(choiceResList)-1])
        if seq_length-1 not in sortChangedList:
            for j in range (len(sortChangedList)):
                aaSequence=aaSequence[:sortChangedList[j]]+str(choiceResList[j])+aaSequence[sortChangedList[j]+1:]
    
    #Substitute with different characteristic residue. 
    #Ex: if the selected residue in Neutral list then the substitute residue from other not Neutral List
    if mutWay == 2:
        choiceResList=[]
        for k in sortChangedList:
            choicePar1=[]
            choicePar2=[]
            if aaSequence[k] in dS.ParPhilic:
                choicePar1 = dS.ParNoPhilic
            if aaSequence[k] in dS.ParPhobic:
                choicePar1 = dS.ParNoPhobic
            if aaSequence[k] in dS.ParNeutral:
                choicePar1 = dS.ParPhilic+dS.ParPhobic
            
            if aaSequence[k] in dS.ParHelix:
                choicePar2 = dS.ParBeta + dS.ParTurn
            if aaSequence[k] in dS.ParBeta:
                choicePar2 = dS.ParHelix + dS.ParTurn
            if aaSequence[k] in dS.ParTurn:
                choicePar2 = dS.ParHelix + dS.ParBeta
            
            choicePar=[]
            for j in range(len(choicePar1)):
                    if choicePar1[j] in choicePar2:
                        choicePar+=[choicePar1[j]]
            if len(choicePar) < 2:
                if aaSequence[k] in dS.ParPhilic:
                    choicePar = dS.ParNoPhilic
                if aaSequence[k] in dS.ParPhobic:
                    choicePar = dS.ParNoPhobic
                if aaSequence[k] in dS.ParNeutral:
                    choicePar = dS.ParPhilic+dS.ParPhobic

            aaWV1=[]
            for j in range(len(choicePar)): 
                for a in range(len(dS.ParAll)):
                    if dS.ParAll[a] == choicePar[j]:
                        aaWV1+=[aaWV[a]]
            choiceRes = "".join(random.choices(population=choicePar, weights=aaWV1, k=1))
            while choiceRes == aaSequence[k]:
                choiceRes = "".join(random.choices(population=choicePar, weights=aaWV1, k=1))
            choiceResList+=[choiceRes]
        
        if seq_length-1 in sortChangedList:
            for j in range(len(sortChangedList)-1):
                aaSequence=aaSequence[:sortChangedList[j]]+str(choiceResList[j])+aaSequence[sortChangedList[j]+1:]
            aaSequence=aaSequence[:(seq_length-1)]+str(choiceResList[len(choiceResList)-1])
        if seq_length-1 not in sortChangedList:
            for j in range (len(sortChangedList)):
                aaSequence=aaSequence[:sortChangedList[j]]+str(choiceResList[j])+aaSequence[sortChangedList[j]+1:]
    return aaSequence

def mutCasesEval(aaSequence, oriSequence, aa_list, seq_length, aaWV):
    print("start eval")
    aaSequence=aaSequence
    aa_list=aa_list
    while aaSequence in aa_list:
        print("Mutated sequence in aa_list")
        aaSequence=oriSequence
        if dI.mutWay == 0:
            mutWay=choice([1, 2])
        aaSequence = mutCases(aaSequence, seq_length, mutWay, mutRes, mutPosition, aaWV)
    aa_list+=[aaSequence]
    valueCV1=dE.NGPSRes(aaSequence)
    valueCV2=dE.RKRes(aaSequence)
    valueCV3=dE.DERes(aaSequence)
    binCV, CanonicalValue=dE.CV(seq_length, valueCV1, valueCV2, valueCV3)
    binII, InstabilityIn=dE.II(seq_length, aaSequence)
    return aaSequence, aa_list, binCV, binII
