#!/usr/bin/python
import os
from random import choice 

#BASED ON IMGT CLASSIFICATION:

##hydrophobic
ParPhobic=["I","V","L","F","C","M","A","W"]
##hydrophilic
ParPhilic=["R","K","E","Q","D","N"]
##Neutral
ParNeutral=["G","T","S","Y","P","H"]
##Positive charge
ParPositive=["R", "H", "K"]
##Negative charge
ParNegative=["D", "E"]

#Linker list BASED ON Argos and George & Heringa
#ParLinker=["T", "S", "P", "G", "D", "K", "Q", "N", "A", "R", "F", "E"]
ParAll=["T", "S", "P", "G", "D", "K", "Q", "N", "A", "V", "E", "R", "I", "Y", "M", "F", "H", "C", "W", "L"]

##Alpha-helix
ParHelix=["M","A","L","E","K"]
##Beta-strand
ParBeta=["Y","F","W","T","V","I"]
##Beta-turn
ParTurn=["G","P","R","D"]
##Polar
ParPolar=["R","N","Q","D","E","H","K","S","T","Y"]
##Non-Polar
ParNonpolar=["A","C","G","I","L","M","F","P","W","V"]


##USED IN THE SCRIPT:
ParCharged=ParPositive + ParNegative
ParNtermin=["Q", "N"]
ParNoQN=list(set(ParAll)-set(ParNtermin))
ParNoQNPhobic=list(set(ParNoQN)-set(ParPhobic))
ParNoQNGPS=list(set(ParAll)-set(["G","P","S"])-set(ParNtermin))
ParNoQNGPSPhobic=list(set(ParNoQNGPS)-set(ParPhobic))
ParNoGPS=list(set(ParAll)-set(["G","P","S"]))
ParNoNGPS=list(set(ParAll)-set(["N","G","P","S"]))
ParNoPhobic=list(set(ParAll)-set(ParPhobic))
ParNoPhilic=list(set(ParAll)-set(ParPhilic))
ParCV=["N","G","P","S","R","K","D","E"]
ParNoCV=list(set(ParAll)-set(ParCV))
ParNoCVPhobic=list(set(ParNoCV)-set(ParPhobic))
ParNoCVPhilPhob=list(set(ParNoCVPhobic)-set(ParPhilic))
ParMultiple=["C","M","W","P","S","V","I","F","Y","L","Q","T","D"]
ParRK=["R","K"]
ParNoGPSPhilic=list(set(ParNoGPS)-set(ParPhilic))
ParNoGPSPhobic=list(set(ParNoGPS)-set(ParPhobic))
ParNeuNoGPS=list(set(ParNeutral)-set(["G","P","S"]))
ParNoGPSPhilPhob=list(set(ParNoGPS)-set(ParPhobic)-set(ParPhilic))
ParNoNGPSPhilPhob=list(set(ParNoNGPS)-set(ParPhobic)-set(ParPhilic))
ParPhilNeu=list(set(ParAll)-set(ParPhobic)-set(ParRK)-set(["N","G","P","S"]))
ParNoNGPSPhobic=list(set(ParNoNGPS)-set(ParPhobic))
ParConservative=["A", "G"]

PhilicHelix = []
PhilicBeta = []
PhilicTurn = []
PhobicHelix = []
PhobicBeta = []
PhobicTurn = []
NeutralHelix = []
NeutralBeta = []
NeutralTurn = []
NegPhilicHelix = []
NegPhilicBeta = []
NegPhilicTurn = []
NegPhobicHelix = []
NegPhobicBeta = []
NegPhobicTurn = []
NegNeutralHelix = []
NegNeutralBeta = []
NegNeutralTurn = []


for k in ParAll:
    if k in ParPhilic and k in ParHelix:
        PhilicHelix+=[k]
    if k in ParPhilic and k in ParBeta:
        PhilicBeta+=[k]
    if k in ParPhilic and k in ParTurn:
        PhilicTurn+=[k]
    
    if k in ParPhobic and k in ParHelix:
        PhobicHelix+=[k]
    if k in ParPhobic and k in ParBeta:
        PhobicBeta+=[k]
    if k in ParPhobic and k in ParTurn:
        PhobicTurn+=[k]
    
    if k in ParNeutral and k in ParHelix:
        NeutralHelix+=[k]
    if k in ParNeutral and k in ParBeta:
        NeutralBeta+=[k]
    if k in ParNeutral and k in ParTurn:
        NeutralTurn+=[k]
    
    if k not in ParPhilic and k in ParHelix:
        NegPhilicHelix+=[k]
    if k not in ParPhilic and k in ParBeta:
        NegPhilicBeta+=[k]
    if k not in ParPhilic and k in ParTurn:
        NegPhilicTurn+=[k]
    
    if k not in ParPhobic and k in ParHelix:
        NegPhobicHelix+=[k]
    if k not in ParPhobic and k in ParBeta:
        NegPhobicBeta+=[k]
    if k not in ParPhobic and k in ParTurn:
        NegPhobicTurn+=[k]
    
    if k not in ParNeutral and k in ParHelix:
        NegNeutralHelix+=[k]
    if k not in ParNeutral and k in ParBeta:
        NegNeutralBeta+=[k]
    if k not in ParNeutral and k in ParTurn:
        NegNeutralTurn+=[k]
    
    

#Function
##all amino acids
def String(length):
    aa_seq=""
    for count in range(length):
        aa_seq+=choice("HRKIFLWAMPVCNGSQYTDE")
    return aa_seq
    
##hydrophobic
def Phobic(length):
    res1=""
    for count in range(length):
        res1+=choice("IVLFCMAW")
    return res1
    
##hydrophilic
def Philic(length):
    res2=""
    for count in range(length):
        res2+=choice("RKEQDN")
    return res2
    
##Neutral
def Neutral(length):
    res3=""
    for count in range(length):
        res3+=choice("GTSYPH")
    return res3
    
##Alpha-helix
def Helix(length):
    res4=""
    for count in range(length):
        res4+=choice("MALEK")
    return res4
    
##Beta-strand
def Beta(length):
    res5=""
    for count in range(length):
        res5+=choice("YFWTVI")
    return res5
    
##Beta-turn
def Turn(length):
    res6=""
    for count in range(length):
        res6+=choice("GPRD")
    return res6

##Polar
def Polar(length):
    res7=""
    for count in range(length):
        res7+=choice("RNQDEHKSTY")
    return res7
    
##Non-Polar
def Nonpolar(length):
    res8=""
    for count in range(length):
        res8+=choice("ACGILMFPWV")
    return res8

##Positive charge
def Positive(length):
    res9=""
    for count in range(length):
        res9+=choice("RHK")
    return res9
    
##Negative charge
def Negative(length):
    res10=""
    for count in range(length):
        res10+=choice("DE")
    return res10

##Charged residue
def Charged(length):
    res11=""
    for count in range(length):
        res11+=choice("RHKDE")
    return res11

##Conservative residue
def Conservative(length):
    res12=""
    for count in range(length):
        res12+=choice("AG")
    return res12
    
##NO in N-terminal residue
def Ntermin(length):
    res13=""
    for count in range(length):
        res13+=choice("HRKIFLWAMPVCGSYTDE")
    return res13
    
