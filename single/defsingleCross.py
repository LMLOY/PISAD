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


def evaluateFindSeq(sequence1, sequence2, seq_length):    
    sameAAList=[]
    for X in range(len(str(sequence1))):
        if sequence1[X] == sequence2[X]:
            sameAAList+=[X]
    print("sameAAList-evaluateFindSeq: "+str(sameAAList))
    redo=""
    if sameAAList == [] or (len(sequence1)-len(sameAAList)) == 1: 
        redo=True
    elif (0 in sameAAList and len(sameAAList) == 1) or (seq_length-1 in sameAAList and len(sameAAList) == 1):
        redo=True
    elif len(sameAAList) == len(str(sequence1)):
        redo="Full"
    else:
        redo=False
    return sequence1, sequence2, sameAAList, redo

def crossProcess (sequence1, sequence2, sameAAList, aaList, seq_length):
    #crossSequence="" 
    sdList=[]
    for X in range(len(sequence1)):
        if X in sameAAList:
            sdList+=["a"]
        else:
            sdList+=["b"]
    print("sdList: "+str(sdList))

    aaCrossList=[]
    for X in range(len(str(sequence1))):
        if sdList[X] == "a":
            aaCrossList+=[sequence1[X]]
        elif X == 0 and sdList[X] == "b":
            selectedOption=choice([sequence1, sequence2])
            aaCrossList+=[selectedOption[X]]        
        elif X!=0 and sdList[X] == "b":
            if sdList[X-1] == "a":
                selectedOption=choice([sequence1, sequence2])
                aaCrossList+=[selectedOption[X]] 
            if sdList[X-1] == "b":
                aaCrossList+=[selectedOption[X]]
        else:
            pass
    crossSequence="".join(aaCrossList)
    print("aaCrossList: "+str(aaCrossList))
    print("crossSequence crossProcess: "+str(crossSequence))
    count = 0
    while (crossSequence in aaList or crossSequence == sequence1 or crossSequence == sequence2) and count < 100:
        aaCrossList=[]
        for X in range(len(str(sequence1))):
            if sdList[X] == "a":
                aaCrossList+=[sequence1[X]]
                #print("1")
            elif X == 0 and sdList[X] == "b":
                selectedOption=choice([sequence1, sequence2])
                aaCrossList+=[selectedOption[X]]
                #print("2")
            elif X!=0 and sdList[X] == "b":
                if sdList[X-1] == "a":
                    selectedOption=choice([sequence1, sequence2])
                    aaCrossList+=[selectedOption[X]] 
                    #print("3")
                if sdList[X-1] == "b":
                    aaCrossList+=[selectedOption[X]]
                    #print("4")
            else:
                pass
        crossSequence="".join(aaCrossList)
        print("aaCrossList: "+str(aaCrossList))
        print("crossSequence crossProcess: "+str(crossSequence))
        count+=1
        print("count-inaaList: "+str(count))
    if count > 99 and crossSequence in aaList:
        print("if99")
        if len(sameAAList) == 1:
            crossSequence1 = sequence1[:sameAAList[0]+1]+sequence2[sameAAList[0]+1:]
            crossSequence2 = sequence2[:sameAAList[0]+1]+sequence1[sameAAList[0]+1:]
            crossSequence=choice([crossSequence1, crossSequence2])
            print("crossSequence crossProcess-if99beforewhile: "+str(crossSequence))  
        else:        
            print("whileif99")
            aaCrossList=[]
            for X in range(len(str(sequence1))):
                if sdList[X] == "a":
                    #selectedOption=sequence1
                    aaCrossList+=[sequence1[X]]
                elif X == 0 and sdList[X] == "b":
                    selectedOption=sequence1
                    aaCrossList+=[selectedOption[X]]        
                elif X!=0 and sdList[X] == "b":
                    if sdList[X-1] == "a":
                        if selectedOption == sequence1:
                            selectedOption=sequence2
                        else:
                            selectedOption=sequence1
                        print("selected: "+str(selectedOption))
                        #selectedOption=choice([sequence1, sequence2])
                        aaCrossList+=[selectedOption[X]] 
                    if sdList[X-1] == "b":
                        aaCrossList+=[selectedOption[X]]
                else:
                    pass   
            print("aaCrossList-while: "+str(aaCrossList))
            crossSequence="".join(aaCrossList)
            print("crossSequence crossProcess-while: "+str(crossSequence))        
    aaList+=[crossSequence]
    return crossSequence, aaList

def crossEval(crosseq, seq_length):
    print("startEval")
    PhilicResidue=dE.PhilicRes(crosseq)
    PhobicResidue=dE.PhobicRes(crosseq)
    valueCV1=dE.NGPSRes(crosseq)
    valueCV2=dE.RKRes(crosseq)
    valueCV3=dE.DERes(crosseq)
    binCV, CanonicalValue=dE.CV(seq_length, valueCV1, valueCV2, valueCV3)
    binII, InstabilityIn=dE.II(seq_length, crosseq)
    return PhilicResidue, PhobicResidue, binCV, CanonicalValue, binII, InstabilityIn