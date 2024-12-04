#!/bin/bash

CODE_DIR="/CODE_DIR" #Directory where all the script is located
PSIPRED_DIR="/PSIPRED_DIR" #Directory where psipred is installed

TYPE=$1 #Initial, Mutation, or Crossover?
POLYPEPTIDE=$2 #True or False? False for single chain only
Start=$3 #Start from. (How many sequences are generated?)
End=$4 #Until. (How many sequences are generated?)

#If running Mutation Process
MUT_DIR="/Users/qiangzhang/Desktop/Jessica实验/single/Mutation" #Input Sequence directory for Mutation process
#If running Crossover Process
CROSS_DIR="/Users/qiangzhang/Desktop/Jessica实验/single/Crossover" #Input Sequence directory for Crossover process

if [[ $# < 4 ]]; then
    echo "Input 3 arguments: Initial/Mutation/Crossover? How many sequences to generate? (in form of integer)"
    echo "Example for generating 50 sequences start from 1: ./GenSeq.sh Initial 1 50"
    echo "or generating 6 mutations per sequence in Mut_Dir: ./GenSeq.sh Mutation 1 6"
    exit 1
fi

if [[ $# > 4 ]]; then
    echo "Too much arguments"
    exit 1
fi

if [[ $1 == "Initial" ]]; then
    echo "Run $1"
    for i in {$Start..$End}; do
        echo "Generating sequence-$i"
        python3 ${CODE_DIR}/$1.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2
        tcsh ${PSIPRED_DIR}/runpsipredplus ${PSIPRED_DIR}/aa_temp.fasta ${PSIPRED_DIR}
        while [[ ! -f ${PSIPRED_DIR}/aa_temp.ss2 ]]; do
            echo "PSIPRED ERROR...Re-Run $1"
            python3 ${CODE_DIR}/appendF.py ${CODE_DIR} #${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2
            python3 ${CODE_DIR}/$1.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2
            tcsh ${PSIPRED_DIR}/runpsipredplus ${PSIPRED_DIR}/aa_temp.fasta ${PSIPRED_DIR}
        done
        python3 ${CODE_DIR}/FilteringSS.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2 ${i} $1
        if [[ $2 == "False" ]]; then
            read -r a b c result < <(tail -n1 ${CODE_DIR}/AllSeqData.txt)
            echo $result
        fi
        if [[ $2 == "True" ]]; then
            read -r a b c d result extra < <(tail -n1 ${CODE_DIR}/AllSeqData.txt)
            echo $result $extra
        fi
        while [[ $result == "F" ]]; do
            echo "Failed. Re-Run $1"
            echo "Generating sequence-" $i
            python3 ${CODE_DIR}/$1.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2
            tcsh ${PSIPRED_DIR}/runpsipredplus ${PSIPRED_DIR}/aa_temp.fasta ${PSIPRED_DIR}
            while [[ ! -f ${PSIPRED_DIR}/aa_temp.ss2 ]]; do
                echo "PSIPRED ERROR...Re-Run $1"
                python3 ${CODE_DIR}/appendF.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2
                python3 ${CODE_DIR}/$1.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2
                tcsh ${PSIPRED_DIR}/runpsipredplus ${PSIPRED_DIR}/aa_temp.fasta ${PSIPRED_DIR}
            done
            python3 ${CODE_DIR}/FilteringSS.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2 ${i} $1
            if [[ $2 == "False" ]]; then
                read -r a b c result < <(tail -n1 ${CODE_DIR}/AllSeqData.txt)
                echo $result
            fi
            if [[ $2 == "True" ]]; then
                read -r a b c d result extra < <(tail -n1 ${CODE_DIR}/AllSeqData.txt)
                echo $result $extra
            fi
        done
    done
fi
if [[ $1 == "Mutation" ]]; then
    for i in ${MUT_DIR}; do
        filename="$(basename "${i%.fasta}")"
        for j in {${Start}..${End}+1}; do 
            python3 ${CODE_DIR}/$1.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2 ${MUT_DIR} ${filename}
            tcsh ${PSIPRED_DIR}/runpsipredplus ${PSIPRED_DIR}/aa_temp.fasta ${PSIPRED_DIR}
            python3 ${CODE_DIR}/FilteringSS.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2 ${filename} $1 ${j}
        done
    done
fi
if [[ $1 == "Crossover" ]]; then
    for i in {${Start}..${End}}; do
        python3 ${CODE_DIR}/$1.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} ${CROSS_DIR} $2
        tcsh ${PSIPRED_DIR}/runpsipredplus ${PSIPRED_DIR}/aa_temp.fasta ${PSIPRED_DIR}
        python3 ${CODE_DIR}/FilteringSS.py ${CODE_DIR} ${CODE_DIR}/Input_Sequence ${PSIPRED_DIR} $2 ${i} $1
    done
fi
