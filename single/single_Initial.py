#!/usr/bin/python
from __future__ import division
import os
import subprocess
import random
from random import choice
import math
import defSeq as dS
import defsingleEquation as dE
import defsingleLogic as dL
import re
import singleParameter as dI

#This code has incorporated with psipred

#Getting Parameter from singleParameter.py
home_dir = dI.home_dir
OutSeq_dir = dI.OutSeq_dir
psipred_dir = dI.psipred_dir
seq_length = dI.seq_length
HHCriteria = dI.HHCriteria
minper_philic = dI.minper_philic
maxper_philic = dI.maxper_philic
maxper_phobic = dI.maxper_phobic
SSCriteria = dI.SSCriteria
Perc_Helix = dI.Perc_Helix
Perc_Beta = dI.Perc_Beta
Perc_Turn = dI.Perc_Turn
aaWV = dI.aaWV


owd = os.getcwd()
os.chdir(OutSeq_dir)

#Inserting a list of previous process sequence
aa_list = dE.aaList()

#Residue Weight Value
if SSCriteria
== "Yes" or SSCriteria == "yes":
	for x
in range(len(dS.ParAll)):
		if dS
	.ParAll[x] in dS.ParHelix:
			aaWV[x] = aaWV[x] + Perc_Helix
				if dS
		.ParAll[x] in dS.ParBeta:
				aaWV[x] = aaWV[x] + Perc_Beta
					if dS
			.ParAll[x] in dS.ParTurn:
					aaWV[x] = aaWV[x] + Perc_Turn
						else
			:
					aaWV[x] = aaWV[x]
						if SSCriteria
				== "No" or SSCriteria == "no":
						aaWV = aaWV

#Normalization of Weight Value
							WVTotal = 0
							for x
					in range(len(aaWV)):
							WVTotal = aaWV[x] + WVTotal
								for x
						in range(len(aaWV)):
								aaWV[x] = aaWV[x] / WVTotal

									print(aaWV)

#Generating sequence
									sequence, aa_list = dL.Logic00(seq_length, aa_list, aaWV)
									print("first: " + str(sequence))

#Count Hydrophilic residue
									sequence, aa_list, binCV, binII = dL.EvaluateSeq(sequence, aa_list, seq_length, aaWV)

##AND Logic for sequence
									while binCV
							== 0 and binII == 0:
									sequence, aa_list = dL.Logic00(seq_length, aa_list, aaWV)
										sequence, aa_list, binCV, binII = dL.EvaluateSeq(sequence, aa_list, seq_length, aaWV)

										while binCV
								== 0 and binII == 1:
										sequence, aa_list = dL.Logic01(aa_list, sequence, seq_length, aaWV)
											sequence, aa_list, binCV, binII = dL.EvaluateSeq(sequence, aa_list, seq_length, aaWV)
											while sequence
									in aa_list and binCV == 0 and binII == 1:
											sequence, aalist = dL.Logic00(seq_length, aa_list, aaWV)
												sequence, aa_list, binCV, binII = dL.EvaluateSeq(sequence, aa_list, seq_length, aaWV)

												while binCV
										== 1 and binII == 0:
												x = 0
													while x
											<1000 + 1 and binCV == 1 and binII == 0:
													sequence = dL.Logic10(sequence)
														sequence, aa_list, binCV, binII = dL.EvaluateSeq(sequence, aa_list, seq_length, aaWV)
														x = x + 1
														print(x)
														if (x > 1000)
												and binCV == 1 and binII == 0:
														sequence, aa_list = dL.Logic00(seq_length, aa_list, aaWV)
															sequence, aa_list, binCV, binII = dL.EvaluateSeq(sequence, aa_list, seq_length, aaWV)

#if binCV == 1 and binII == 1:
#os.chdir(psipred_dir)
#f= open("aa_temp.fasta","w+")
#f.write(">aaTemp \n")
#f.write("".join(sequence))
#os.chdir(owd)
#os.chdir(home_dir)
#with open("AllSeqData.txt", "a") as file:
