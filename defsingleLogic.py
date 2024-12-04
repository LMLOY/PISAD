#!/usr/bin/python
from __future__ import division
import os
import random
from random import choice
import defSeq as dS
import defsingleEquation as dE
import re
import singleParameter as dI

minper_philic=dI.minper_philic
maxper_philic=dI.maxper_philic
maxper_phobic=dI.maxper_phobic
HHCriteria = dI.HHCriteria

#Criteria
##25% < Hydrophilic aa (min: 12) < 75%, Hydrophobic aa (max: 8) <50%
##at least one charged aa for every 5 aa. Start from residue-4 (4, 9, 14, 19, 24)
##no Q and N at first four residues
##no combination of D with G, P, and S

constantA = 0.03
constantB = 15.43
constantC = 29.56
constantD = 1.71

def Logic00(seq_length, aa_list, aaWV):
    #Generating sequence
    aa_list=aa_list
    sequence=dE.seq(seq_length, aaWV)
    while sequence in aa_list:
        print("same sequence")
        sequence=dE.seq(seq_length, aaWV)
    aa_list+=[sequence]
    return sequence, aa_list

def EvaluateSeq(sequence, aa_list, seq_length, aaWV):
    aa_list = aa_list
    if HHCriteria == "Yes" or HHCriteria == "yes":
        #Evaluate Hydrophilic residue
        PhilicResidue = dE.PhilicRes(sequence)
        if PhilicResidue < int(minper_philic * seq_length):
            sequence = dE.changePhilicMin(seq_length, sequence, aaWV)
        elif PhilicResidue > int(maxper_philic * seq_length):
            sequence = dE.changePhilicMax(seq_length, sequence, aaWV)
        else:
            sequence = sequence
        if sequence in aa_list:
            pass
        else:         
            aa_list+=[sequence]
        
        #Evaluate Hydrophobic residue
        PhobicResidue = dE.PhobicRes(sequence)
        if PhobicResidue > int(maxper_phobic * seq_length):
            sequence=dE.changePhobic(sequence, seq_length, aaWV)
        else:
            sequence=sequence
        
        #Final Sequence
        PhilicResidue = dE.PhilicRes(sequence)
        PhobicResidue = dE.PhobicRes(sequence)
        if sequence in aa_list:
            pass
        else:         
            aa_list+=[sequence]
    if HHCriteria == "No" or HHCriteria == "no":
        sequence = sequence
    
    ##Calculating canonical value (Davis et al. 1999)
    valueCV1=dE.NGPSRes(sequence)
    valueCV2=dE.RKRes(sequence)
    valueCV3=dE.DERes(sequence)
    print("CV1: "+str(valueCV1)+" CV2: " +str(valueCV2)+ " CV3: " +str(valueCV3))
    binCV, CanonicalValue = dE.CV(seq_length, valueCV1, valueCV2, valueCV3)

    ##Calculating instability index (Guruparasad et al. 1990)
    binII, InstabilityIn=dE.II(seq_length, sequence)
    return sequence, aa_list, binCV, binII

def minDiffDE(valueCV1, length): #depend on DE
    fixA= constantA * length
    fixB= constantB * valueCV1 / constantC
    fixC= constantD * length / constantC
    nCV1 = int(length / 2)
    nLeft = length - nCV1
    minCV2CV3=False
    for jCV2 in range (nLeft):
        for jCV3 in range (nLeft):
            if minCV2CV3 == False:
                if jCV2 + jCV3 < nLeft+1:
                    mindiff = jCV2-jCV3
                    valueA = (-1) * abs(mindiff - fixA)
                    valueB = fixC - fixB
                    if valueA < valueB:
                        minCV2CV3=True
    print("mindiff: "+str(mindiff))
    #print(minCV2CV3)
    return mindiff

def minDiffRK(valueCV1, length):
    fixA= constantA * length
    fixB= constantB * valueCV1 / constantC
    fixC= constantD * length / constantC
    nCV1 = int(length / 2)
    nLeft = length - nCV1
    minCV2CV3=False
    for jCV2 in range (nLeft):
        for jCV3 in range (nLeft):
            if minCV2CV3 == False:
                if jCV2 + jCV3 < nLeft+1:
                    mindiff = jCV3-jCV2
                    valueA = (-1) * abs(mindiff - fixA)
                    valueB = fixC - fixB
                    if valueA < valueB:
                        minCV2CV3=True
    print("mindiffRK: "+str(mindiff))
    #print(minCV2CV3)
    return mindiff

def recall(seq_length, sequence, oriseq, parameter, aaWV):
    aaPhob=sequence
    valueCV1=dE.NGPSRes(aaPhob) 
    valueCV2=dE.RKRes(aaPhob)
    valueCV3=dE.DERes(aaPhob)
    mindiff=minDiffDE(valueCV1, seq_length)
    #mindiffRK=minDiffRK(valueCV1)
    diffCV3n2=valueCV3-valueCV2
    diffCV2n3=valueCV2-valueCV3
    PhilicResidue=dE.PhilicRes(aaPhob)
    diff=abs(diffCV3n2 - abs(mindiff))
    totalphilic=diff+PhilicResidue
    condition = ""
    
    if diff%2 == 0 and diffCV2n3 > int(diff/2)-1:
        XChanged_List=[]
        for X in range(len(str(aaPhob))):
            if aaPhob[X] in dS.ParRK:
                XChanged_List+=[X]
        aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, int(diff/2), XChanged_List, dS.ParNegative, aaWV)
    elif diff%2 != 0 and diffCV2n3 > (int(diff/2)+diff%2)-1:
        XChanged_List=[]
        for X in range(len(str(aaPhob))):
            if aaPhob[X] in dS.ParRK:
                XChanged_List+=[X]
        aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, int(diff/2), XChanged_List, dS.ParNegative, aaWV)
        XChanged_List=[]
        for X in range(len(str(aaPhob))):
            if aaPhob[X] in dS.ParRK and X not in indexAA:
                XChanged_List+=[X]
        aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diff%2, XChanged_List, dS.ParPhilNeu, aaWV)
    else:
        halfNum=PhilicResidue-int(diff/2)-diff%2
        if HHCriteria == "Yes" or HHCriteria == "yes":
            condition = halfNum > int(minper_philic*seq_length)-1 and (valueCV2 > (int(diff/2)+diff%2)-1) and (valueCV2>valueCV3)
        if HHCriteria == "No" or HHCriteria == "no":
            condition = (valueCV2 > (int(diff/2)+diff%2)-1) and (valueCV2>valueCV3)
        if condition:
            RKChanged_List=[]
            for X in range(len(str(aaPhob))):
                if aaPhob[X] in dS.ParRK:
                    RKChanged_List+=[X]
            if len(RKChanged_List) > diff-1:
                aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diff, RKChanged_List, dS.ParPhilNeu, aaWV)
            else:
                if len(RKChanged_List) - (int(diff/2)+diff%2) > -1:
                    aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, int(diff/2), RKChanged_List, dS.ParNegative, aaWV)
                    if diff%2 != 0:
                        XChanged_List=[]
                        for X in range(len(str(aaPhob))):
                            if X in RKChanged_List and X not in indexAA:
                                XChanged_List+=[X]
                        aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diff%2, XChanged_List, dS.ParPhilNeu, aaWV)
                else:
                    XChanged_List=[]
                    for X in range(len(str(aaPhob))):
                        if aaPhob[X] in dS.ParRK:
                            XChanged_List+=[X]
                    print(XChanged_List)
                    aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, int(diff/2), XChanged_List, dS.ParNegative, aaWV)
                    valueCV2=dE.RKRes(aaPhob)
                    valueCV3=dE.DERes(aaPhob)
                    different=valueCV3-valueCV2
                    needtochange=abs(mindiff-different)
                    if needtochange != 0:
                        if valueCV2 > needtochange:
                            XChanged_List=[]
                            for X in range(len(str(aaPhob))):
                                if aaPhob[X] in dS.ParRK:
                                    XChanged_List+=[X]
                            aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, needtochange, XChanged_List, dS.ParPhilNeu, aaWV)
                        else:
                            if valueCV2 > int(needtochange/2)+needtochange%2:
                                XChanged_List=[]
                                for X in range(len(str(aaPhob))):
                                    if aaPhob[X] in dS.ParRK:
                                        XChanged_List+=[X]
                                aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, int(needtochange/2)+needtochange%2, XChanged_List, dS.ParNegative, aaWV)
                            else:
                                aaPhob =oriseq
        elif totalphilic < parameter:
            XChanged_List=[]
            for X in range(len(str(aaPhob))):
                if aaPhob[X] not in dS.ParCV:
                    XChanged_List+=[X]
            aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diff, XChanged_List, dS.ParNegative, aaWV)
        else:
            aaPhob=oriseq
    return aaPhob
    
def Logic01(aa_list, sequence, seq_length, aaWV):
    aaPhob = sequence
    oriseq = sequence
    valueCV1 = dE.NGPSRes(aaPhob) 
    valueCV2 = dE.RKRes(aaPhob)
    valueCV3 = dE.DERes(aaPhob)
    mindiff = minDiffDE(valueCV1, seq_length)
    #mindiffRK=minDiffRK(valueCV1)
    diffCV3n2 = valueCV3-valueCV2
    diffCV2n3 = valueCV2-valueCV3
    PhilicResidue = dE.PhilicRes(aaPhob)
    diff = abs(diffCV3n2 - abs(mindiff))
    totalphilic = diff+PhilicResidue
    #need=abs(diffCV3n2 + abs(mindiff)) 
    if HHCriteria == "Yes" or HHCriteria == "yes":
        if maxper_philic != []:
            parameter=int(maxper_philic * seq_length)+1
            aaPhob=recall(seq_length, aaPhob, oriseq, parameter, aaWV)
            aaPhob, aa_list, binCV, binII, = EvaluateSeq(aaPhob, aa_list, seq_length, aaWV)    
            
            if binCV == 0 and binII == 1:       
                if valueCV1 > int(0.2*seq_length):
                    diffCV1=valueCV1-int(0.2*seq_length)
                    XChanged_List=[]
                    for X in range(len(str(aaPhob))):
                       if aaPhob[X] in ["N","G","P","S"]:
                           XChanged_List+=[X]
                    aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diffCV1, XChanged_List, list(set(dS.ParNoNGPS)-set(["Q"])), aaWV)
                    PhobicResidue=dE.PhobicRes(aaPhob)
                    if PhobicResidue > int(maxper_phobic*seq_length):
                        diffPhob=PhobicResidue-int(maxper_phobic*seq_length)
                        XChanged_List=[]
                        for X in range(len(str(aaPhob))):
                            if X in indexAA and aaPhob[X] in dS.ParPhobic:
                                XChanged_List+=[X]
                        aaPhob, indexAA1=dE.changingseq(aaPhob, seq_length, diffPhob, XChanged_List, dS.ParNoNGPSPhobic, aaWV)
                        indexAA=indexAA+indexAA1
                    PhilicResidue=dE.PhilicRes(aaPhob)
                    if PhilicResidue > int(maxper_philic*seq_length):
                        diffPhil=PhilicResidue-int(maxper_philic*seq_length)
                        XChanged_List=[]
                        for X in range(len(str(aaPhob))):
                            if X in indexAA and aaPhob[X] in dS.ParPhilic:
                                XChanged_List+=[X]
                        aaPhob, indexAA2=dE.changingseq(aaPhob, seq_length, diffPhil, XChanged_List, dS.ParNoNGPSPhilPhob, aaWV)
                        indexAA=indexAA+indexAA2
                    aaPhob=recall(seq_length, aaPhob, oriseq, parameter, aaWV)
                else: 
                    changeCV1=valueCV1
                    while mindiff != 0:
                        changeCV1=changeCV1
                        mindiff=minDiffDE(changeCV1, length)
                        changeCV1=changeCV1-1
                    diffCV1=valueCV1-changeCV1
                    if diffCV1 > -1:
                        XChanged_List=[]
                        for X in range(len(str(aaPhob))):
                           if aaPhob[X] in ["N","G","P","S"]:
                               XChanged_List+=[X]
                        aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diffCV1, XChanged_List, list(set(dS.ParNoNGPS)-set(["Q"])), aaWV)
                        PhobicResidue=dE.PhobicRes(aaPhob)
                        if PhobicResidue > int(maxper_phobic*seq_length):
                            diffPhob=PhobicResidue-int(maxper_phobic*seq_length)
                            XChanged_List=[]
                            for X in range(len(str(aaPhob))):
                                if X in indexAA and aaPhob[X] in dS.ParPhobic:
                                    XChanged_List+=[X]
                            aaPhob, indexAA1=dE.changingseq(aaPhob, seq_length, diffPhob, XChanged_List, dS.ParNoNGPSPhobic, aaWV)
                            indexAA=indexAA+indexAA1
                        PhilicResidue=dE.PhilicRes(aaPhob)
                        if PhilicResidue > int(maxper_philic*seq_length):
                            diffPhil=PhilicResidue-int(maxper_philic*seq_length)
                            XChanged_List=[]
                            for X in range(len(str(aaPhob))):
                                if X in indexAA and aaPhob[X] in dS.ParPhilic:
                                    XChanged_List+=[X]
                            aaPhob, indexAA2=dE.changingseq(aaPhob, seq_length, diffPhil, XChanged_List, dS.ParNoNGPSPhilPhob, aaWV)
                            indexAA=indexAA+indexAA2
                        aaPhob=recall(seq_length, aaPhob, oriseq, parameter, aaWV)
                    else:
                        aaPhob=oriseq
                        aa_list=aa_list
        else:
            parameter=seq_length-valueCV1
            aaPhob=recall(sesq_length, aaPhob, oriseq, parameter, aaWV)
            aaPhob, aa_list, binCV, binII = EvaluateSeq(aaPhob, aa_list, seq_length, aaWV)    
            
            if binCV == 0 and binII == 1:   
                if valueCV1 > int(0.2*seq_length):
                    diffCV1=valueCV1-int(0.2*seq_length)
                    XChanged_List=[]
                    for X in range(len(str(aaPhob))):
                       if aaPhob[X] in ["N","G","P","S"]:
                           XChanged_List+=[X]
                    aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diffCV1, XChanged_List, list(set(dS.ParNoNGPS)-set(["Q"])), aaWV)
                    PhobicResidue=dE.PhobicRes(aaPhob)
                    if PhobicResidue > int(maxper_phobic*seq_length):
                        diffPhob=PhobicResidue-int(maxper_phobic*seq_length)
                        XChanged_List=[]
                        for X in range(len(str(aaPhob))):
                            if X in indexAA and aaPhob[X] in dS.ParPhobic:
                                XChanged_List+=[X]
                        aaPhob, indexAA1=dE.changingseq(aaPhob, seq_length, diffPhob, XChanged_List, dS.ParNoNGPSPhobic, aaWV)
                        indexAA=indexAA+indexAA1
                    aaPhob=recall(seq_length, aaPhob, oriseq, parameter, aaWV)
                    
                else: 
                    changeCV1=valueCV1
                    while mindiff != 0:
                        changeCV1=changeCV1
                        mindiff=minDiffDE(changeCV1, seq_length)
                        changeCV1=changeCV-1
                    diffCV1=valueCV1-changeCV1
                    if diffCV1 > -1:
                        XChanged_List=[]
                        for X in range(len(str(aaPhob))):
                           if aaPhob[X] in ["N","G","P","S"]:
                               XChanged_List+=[X]
                        aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diffCV1, XChanged_List, list(set(dS.ParNoNGPS)-set(["Q"])), aaWV)
                        PhobicResidue=dE.PhobicRes(aaPhob)
                        if PhobicResidue > int(maxper_phobic*seq_length):
                            diffPhob=PhobicResidue-int(maxper_phobic*seq_length)
                            XChanged_List=[]
                            for X in range(len(str(aaPhob))):
                                if X in indexAA and aaPhob[X] in dS.ParPhobic:
                                    XChanged_List+=[X]
                            aaPhob, indexAA1=dE.changingseq(aaPhob, seq_length, diffPhob, XChanged_List, dS.ParNoNGPSPhobic, aaWV)
                            indexAA=indexAA+indexAA1
                        aaPhob=recall(seq_length, aaPhob, oriseq, parameter, aaWV)
                    else:
                        aaPhob=oriseq
                        aa_list=aa_list
    if HHCriteria == "No" or HHCriteria == "no":
        parameter=seq_length-valueCV1
        aaPhob=recall(sesq_length, aaPhob, oriseq, parameter, aaWV)
        aaPhob, aa_list, binCV, binII = EvaluateSeq(aaPhob, aa_list, seq_length, aaWV)    
        if binCV == 0 and binII == 1:   
            if valueCV1 > int(0.2*seq_length):
                diffCV1=valueCV1-int(0.2*seq_length)
                XChanged_List=[]
                for X in range(len(str(aaPhob))):
                   if aaPhob[X] in ["N","G","P","S"]:
                       XChanged_List+=[X]
                aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diffCV1, XChanged_List, list(set(dS.ParNoNGPS)-set(["Q"])), aaWV)
                aaPhob=recall(seq_length, aaPhob, oriseq, parameter, aaWV)
            else: 
                changeCV1=valueCV1
                while mindiff != 0:
                    changeCV1=changeCV1
                    mindiff=minDiffDE(changeCV1, seq_length)
                    changeCV1=changeCV-1
                diffCV1=valueCV1-changeCV1
                if diffCV1 > -1:
                    XChanged_List=[]
                    for X in range(len(str(aaPhob))):
                       if aaPhob[X] in ["N","G","P","S"]:
                           XChanged_List+=[X]
                    aaPhob, indexAA=dE.changingseq(aaPhob, seq_length, diffCV1, XChanged_List, list(set(dS.ParNoNGPS)-set(["Q"])), aaWV)
                    aaPhob=recall(seq_length, aaPhob, oriseq, parameter, aaWV)
                else:
                    aaPhob=oriseq
                    aa_list=aa_list
    return aaPhob, aa_list

def countDuplicate(sequence):
    TLet=dE.duplicate(sequence, "T")
    SLet=dE.duplicate(sequence, "S")
    PLet=dE.duplicate(sequence, "P")
    GLet=dE.duplicate(sequence, "G")
    DLet=dE.duplicate(sequence, "D")
    KLet=dE.duplicate(sequence, "K")
    QLet=dE.duplicate(sequence, "Q")
    NLet=dE.duplicate(sequence, "N")
    ALet=dE.duplicate(sequence, "A")
    RLet=dE.duplicate(sequence, "R")
    FLet=dE.duplicate(sequence, "F")
    ELet=dE.duplicate(sequence, "E")
    HLet=dE.duplicate(sequence, "H")
    ILet=dE.duplicate(sequence, "I")
    LLet=dE.duplicate(sequence, "L")
    WLet=dE.duplicate(sequence, "W")
    MLet=dE.duplicate(sequence, "M")
    VLet=dE.duplicate(sequence, "V")
    CLet=dE.duplicate(sequence, "C")
    YLet=dE.duplicate(sequence, "Y")
    
    countTLet=len(TLet)
    countSLet=len(SLet)
    countPLet=len(PLet)
    countGLet=len(GLet)
    countDLet=len(DLet)
    countKLet=len(KLet)
    countQLet=len(QLet)
    countNLet=len(NLet)
    countALet=len(ALet)
    countRLet=len(RLet)
    countFLet=len(FLet)
    countELet=len(ELet)                  
    countHLet=len(HLet)
    countILet=len(ILet)
    countLLet=len(LLet)
    countWLet=len(WLet)
    countMLet=len(MLet)
    countVLet=len(VLet)
    countCLet=len(CLet)
    countYLet=len(YLet)
         
    return countRLet, countKLet, countDLet, countELet,countHLet, countILet, countFLet, countLLet, countWLet, countALet, countMLet, countPLet, countVLet, countCLet, countNLet, countGLet, countSLet, countQLet, countYLet, countTLet

def factorial(n):
    if n == 0:
        fact=1
    else:
        fact=1
        for i in range (1, n+1):
            fact=i*fact
    return fact
    
def letterCalculate(sequence, length):
    countRLet, countKLet, countDLet, countELet,countHLet, countILet, countFLet, countLLet, countWLet, countALet, countMLet, countPLet, countVLet, countCLet, countNLet, countGLet, countSLet, countQLet, countYLet, countTLet=countDuplicate(sequence)
    print(countRLet, countKLet, countDLet, countELet,countHLet, countILet, countFLet, countLLet, countWLet, countALet, countMLet, countPLet, countVLet, countCLet, countNLet, countGLet, countSLet, countQLet, countYLet, countTLet)
    seqfact=factorial(length)
    Tfact=factorial(countTLet)
    Sfact=factorial(countSLet)
    Pfact=factorial(countPLet)
    Gfact=factorial(countGLet)
    Dfact=factorial(countDLet)
    Kfact=factorial(countKLet)
    Qfact=factorial(countQLet)
    Nfact=factorial(countNLet)
    Afact=factorial(countALet)
    Rfact=factorial(countRLet)
    Ffact=factorial(countFLet)
    Efact=factorial(countELet)
    Hfact=factorial(countHLet)
    Ifact=factorial(countILet)
    Lfact=factorial(countLLet)
    Wfact=factorial(countWLet)
    Mfact=factorial(countMLet)
    Vfact=factorial(countVLet)
    Cfact=factorial(countCLet)
    Yfact=factorial(countYLet)
    totalfact=seqfact/(Rfact*Kfact*Dfact*Efact*Hfact*Ifact*Ffact*Lfact*Wfact*Afact*Mfact*Pfact*Vfact*Cfact*Nfact*Gfact*Sfact*Qfact*Yfact*Tfact)
    print(totalfact)
    return totalfact
   
def Logic10(sequence):
    shuffleList=[]
    for Z in range(len(str(sequence))):
        shuffleRes=sequence[Z]
        shuffleList+=[shuffleRes]
    aaShuffle="".join(random.sample(shuffleList,len(shuffleList)))
    return aaShuffle
