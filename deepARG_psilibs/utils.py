#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio import SeqIO
import numpy as np
import os
import fileinput
import re

def toFile(seq, seq_id, out_dir):
    file_path = os.path.join("{}".format(out_dir),
                             "{}.fasta".format(seq_id.split("|")[0]))
    outfile = open(file_path, "w+")
    header = ">" + seq_id + "\n"
    outfile.write(header)
    outfile.write(str(seq) + "\n")
    outfile.close()
    return str(file_path)

def psiblast(query, db, out_dir, out_file):
    out_file = os.path.join(out_dir, "{}.pssm".format(out_file))
    psi_cline = NcbipsiblastCommandline('psiblast', db = db,
                                        query = query, evalue = 0.001,
                                        num_iterations = 3,
                                        save_pssm_after_last_round=True,
                                        out_ascii_pssm = out_file )

    psi_cline()
    return out_file

def readToMatrix(input_matrix):
    #print "start to read PSSM matrix"
    PSSM = []
    p = re.compile(r'-*[0-9]+')
    for line, strin in enumerate(input_matrix):
        if line > 2:
            str_vec = []
            overall_vec = strin.split()
            if len(overall_vec) == 0:
                break
            str_vec.extend(overall_vec[1])
            if(len(overall_vec) < 44):
                print("There is a mistake in the pssm file")
                print("Try to correct it")
                for cur_str in overall_vec[2:]:
                    str_vec.extend(p.findall(cur_str))
                    if(len(str_vec) >= 21):
                        if(len(str_vec)) >21:
                            exit(1)
                        break;
                print("Done")
            else:
                str_vec = strin.split()[1:42]
            if len(str_vec) == 0:
                break
            #str_vec_positive=map(int, str_vec[1:])
            PSSM.append(str_vec)
    fileinput.close()
    #print "finish to read PSSM matrix"
    PSSM = np.array(PSSM)
    return PSSM


def average(matrixSum, seqLen):
    # average the summary of rows
    matrix_array = np.array(matrixSum)
    matrix_array = np.divide(matrix_array, seqLen)
    matrix_array_shp = np.shape(matrix_array)
    matrix_average = [(np.reshape(matrix_array, (matrix_array_shp[0] * matrix_array_shp[1], )))]
    return np.array(matrix_average)

def handleRows(PSSM, SWITCH, COUNT):
    '''
    if SWITCH=0, we filter no element.
    if SWITCH=1, we filter all the negative elements.
    if SWITCH=2, we filter all the negative and positive elements greater than expected.
    '''
    '''
    if COUNT=20, we generate a 20-dimension vector.
    if COUNT=400, we generate a 400-dimension vector.
    '''
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"

    matrix_final = [ [0.0] * 20 ] * int(COUNT/20)
    matrix_final=np.array(matrix_final)
    seq_cn = 0

    PSSM_shape=np.shape(PSSM)
    for i in range(PSSM_shape[0]):
        seq_cn += 1
        str_vec=PSSM[i]
        str_vec_positive=list(map(int, str_vec[1:21]))
        str_vec_positive=np.array(str_vec_positive)
        if SWITCH==1:
            str_vec_positive[str_vec_positive<0]=0
        elif SWITCH==2:
            str_vec_positive[str_vec_positive<0]=0
            str_vec_positive[str_vec_positive>7]=0
        if COUNT==20:
            matrix_final[0]=list(map(sum, list(zip(str_vec_positive, matrix_final[0]))))
        elif COUNT==400:
            if Amino_vec.find(str_vec[0]) > 0:
                matrix_final[Amino_vec.index(str_vec[0])] = list(map(sum, list(zip(str_vec_positive, matrix_final[Amino_vec.index(str_vec[0])]))))

    return matrix_final

