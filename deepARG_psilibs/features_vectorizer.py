#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

from sklearn.model_selection import train_test_split
from PyBioMed.PyProtein import CTD
from Bio import SeqIO
import numpy as np
import fileinput
from deepARG_psilibs.utils import *
import pickle
import os
import csv

class Features_extractor():
    def __init__(self,
                 protein_fasta_file,
                 project_root):
        self.protein_sequence = protein_fasta_file
        self.pr = project_root
        self.amr_gene_class_dict = {
           "aminoglycoside": 1, "acriflavin": 2, "acriflavine": 3,
           "aminocoumarin": 4, "bacitracin": 5, "beta_lactam": 6,
           "bleomycin": 7, "chloramphenicol": 8, "deoxycholate": 9,
           "doxorubicin": 10, "multidrug": 11, "elfamycin": 12,
           "ethambutol": 13, "fluoramphenicol": 14, "fosfomycin": 15,
           "fosmidomycin": 16, "fusaric_acid": 17, "fusidic_acid": 18,
           "quinolone": 19, "glycopeptide": 20, "tetracycline": 21,
           "isoniazid": 22, "kasugamycin": 23, "linezolid": 24,
           "lipopeptide": 25, "mupirocin": 26, "nitrofuratoin": 27,
           "peptide": 28, "phenicol": 29, "pleuromutilin": 30,
           "puromycin": 31, "tunicamycin": 32, "macrolide-lincosamide-streptogramin": 33,
           "polymyxin": 34, "pyrazinamide": 35, "rifampin": 36,
           "roxithromycin": 37, "t_chloride": 38,
           "qa_compound": 40, "tetracenomycin": 41, "viomycin": 42,
           "streptothricin": 43, "sulfonamide": 44, "polyamine": 45

           }

    def label_AMR_gene(self, fasta_header):
        entries = fasta_header.split("|")
        amr_class = entries[3]
        return self.amr_gene_class_dict.get(amr_class, 0)
    
    def CTD_feature_extractor(self, protein_sequence):
        value_list = []
        features_list = ["_Polarizability", "_SolventAccessibility", "_SecondaryStr", 
                       "_Charge", "_Polarity", "_NormalizedVDWV", "_Hydrophobicity"]
        ctd_features = CTD.CalculateCTD(protein_sequence)
        for feature_keys in ctd_features.keys():
            value = ctd_features.get(feature_keys)
            if value > 1:
                value = round((value/100), 3)
            value_list.append(value)
        return np.array(value_list, dtype=np.float32)
    
    def get_data(self):
        target_list = []
        initial_array = ""

        for record in SeqIO.parse(self.protein_sequence, "fasta"):
            target_list.append(self.label_AMR_gene(record.id))
            if initial_array == "":
                initial_array = np.array(self.CTD_feature_extractor(str(record.seq)), dtype=np.float32)
            else:
                feature_array = self.CTD_feature_extractor(str(record.seq))
                initial_array = np.vstack((initial_array, feature_array))
        target_array = np.array(target_list, dtype=np.int32)
        feature_tuple = (initial_array, target_array)

        return feature_tuple
    
    def PSSM_comp_Extractor(self):
        target_list = []
        initial_array = ""

        for record in SeqIO.parse(self.protein_sequence, "fasta"):
            file_path = toFile(str(record.seq),
                               record.id,
                               "{}/temp".format(self.pr))
            pssm_path = psiblast(file_path,
                                 "{}/database/card_ardb.fasta".format(self.pr),
                                 "{}/pssm_files/".format(self.pr),
                                 record.id.split("|")[0])
            if os.path.exists(pssm_path):
                pssm_array = readToMatrix(fileinput.input(pssm_path))
                final_array = average(handleRows(pssm_array, 0, 400),
                                  len(pssm_array))
                target_list.append(self.label_AMR_gene(record.id))
                if initial_array == "":
                    initial_array = final_array
                else:
                    initial_array = np.vstack((initial_array, final_array))
                os.remove(file_path)
                os.remove(pssm_path)
            else:
                print("{} was not created".format(pssm_path))
        target_array = np.array(target_list, dtype=np.int32)
        feature_tuple = (initial_array, target_array)

        return feature_tuple
    
    def PSSM_comp_Extractor_predict(self, model_path, out):
        loaded_model = pickle.load(open(model_path, 'rb'))
        with open(out, mode='w+') as out_file:
            prediction_writer = csv.writer(out_file,
                                           delimiter=',',
                                           quotechar='"',
                                           quoting=csv.QUOTE_MINIMAL)
            prediction_writer.writerow(["Sequence ID",
                                        "Sequence Length",
                                        "Predicted ABR category",
                                        "Probability"])
            for record in SeqIO.parse(self.protein_sequence, "fasta"):
                file_path = toFile(str(record.seq),
                                   record.id,
                                   "{}/temp/".format(self.pr))
                pssm_path = psiblast(file_path,
                                     "{}/database/card_ardb.fasta".format(self.pr),
                                     "{}/pssm_files/".format(self.pr),
                                     record.id.split("|")[0])
                if os.path.exists(pssm_path):
                    pssm_array = readToMatrix(fileinput.input(pssm_path))
                    final_array = average(handleRows(pssm_array, 0, 400),
                                      len(pssm_array))
                    prediction = int("".join([str(i) for i in loaded_model.predict(final_array).tolist()]))
                    pred_prob = max(loaded_model.predict_proba(final_array)[0])
                    pred_category = list(self.amr_gene_class_dict.keys())[list(self.amr_gene_class_dict.values()).index(prediction)]
                    entry = [record.description,
                             str(len(record.seq)),
                             pred_category,
                             pred_prob]
                    prediction_writer.writerow(entry)




