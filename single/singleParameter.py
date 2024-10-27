#!/usr/bin/python
from __future__ import division 
import os 
import random 
import sys 
from random import choice 
import defSeq as dS 

#Directory
#home_dir = '/home/dinglab/Jessica' #Parent directory where all the codes are located
#OutSeq_dir = '/home/dinglab/Jessica/Sequence' #output sequence directory
#psipred_dir = '/home/dinglab/psipred' #PSIPRED installed directory
#Mut_dir = '/home/dinglab/Jessica/Mutation' #Input sequence for mutation
#Cross_dir = '/home/dinglab/Jessica/Crossover' #Input sequence for crossover

home_dir = '/Users/qiangzhang/Desktop/Jessica实验/single' 
#Parent directory where all the codes are located
OutSeq_dir = 'Sequence' 
#output sequence directory
psipred_dir = '/Users/qiangzhang/Desktop/Jessica实验/single' 
#PSIPRED installed directory
Mut_dir = 'Mutation' 
#Input sequence for mutation
Cross_dir = 'Crossover' 
#Input sequence for crossover

seq_length = 20 
#sequence length

#Input for filtering the sequence based on secondary structure
coilPercentage = 0.3 
#percentage of total coil in a sequence 

#Set Hydrophobic/Hdyrophilic Residue
##25% < Hydrophilic aa < 75%, Hydrophobic aa < 50%
HHCriteria = "Yes" 
#want to set or not? Yes or No.
minper_philic = 0.25 
#minimum percentage of hydrophilic residue
maxper_philic = 0.50 
#maximum percentage of hydrophilic residue (optional). If don't want to specify the maximum percentage of hydrophilic residue, input 0
maxper_phobic = 0.35 
#maximum percentage of hydrophobic residue

#Set Secondary Structure Residue Criteria
SSCriteria = "Yes" 
#want to set or not? Yes or No.
Perc_Helix = 0.45 
#increase weight value to choose amino acid that commonly construct helix
Perc_Beta = 0.35 
#increase weight value to choose amino acid that commonly construct beta
Perc_Turn = 0.15 
#increase weight value to choose amino acid that commonly construct turn
#the weight value is increased to the aaWV in the algorithm when run.

#Set Weight Value (probablility) for selection of list of amino acid below.
#["T", "S", "P", "G", "D", "K", "Q", "N", "A", "V", "E", "R", "I", "Y", "M", "F", "H", "C", "W", "L"]
aaWV =[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] 
#if no amino acid preferences, set all to 1


#Input for mutated sequence
Consecutive = "Yes" 
#0 for not sequential position mutation, 1 for sequential, and 2 for randomly choose between them.
mutPosition =[] 
#if mutHow=1 or 2, you can specify the first mutation position. if not, then input [].
mutWay = 0 
#Which way to do mutation (optional).  If don't want to specify it, input 0.
#Currently there are 5 ways. 
#If mutHow = 0, then the choices are 1, 2, 3. 
#If mutHow = 1. then the choices are 4, 5.
mutRes = 3 
#Number of residue/s mutated in a sequence


