# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 13:30:20 2022

@author: user
helper script subsetting 5` and 3` UTR regions in specific way

original file:
# Column 1: chromosome name
# Column 2: source (First author and publication year, assigned by SGD)
# Column 3: feature type (Sequence ontology term assigned by SGD, where applicable)
# Column 4: start coordinate of feature (1-based)
# Column 5: end coordinate of feature (end inclusive)
# Column 6: score where applicable, "." otherwise
# Column 7: strand
# Column 8: frame where applicable, "." otherwise
# Column 9: Feature attributes, detailed below,

final file:
# Column 1: gene name
# Column 2: chromosome name
# Column 3: start position
# Column 4: end position
"""
import csv
path_to_comp_file =  "C:/Users/user/Desktop/BPC_FP/programming/5_UTR/Nagalakshmi_2008_3UTRs_V64.gff3"
path_to_easy_file =  "C:/Users/user/Desktop/BPC_FP/programming/5_UTR/3_UTR.gff3"

def get_ID_name(string):
    features = string.split(';')
    ID = features[0].split('=')
    return ID[1]
    

def gff_easy_format(in_path, out_path):
    with open(out_path, 'w', newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')    
        with open(in_path) as f:
            for line in f:
                if not line.lstrip().startswith('#'):
                    words = line.split()
                    if len(words) != 0:
                        new_line = [get_ID_name(words[8]), words[0], words[3], words[4]]
                        tsv_writer.writerow(new_line)

#gff_easy_format(path_to_comp_file, path_to_easy_file)


"""
Created on Mon Jul 06 13:30:20 2022

@author: user
get features of a transcript in a specific way, so the data can be used as the 
in_features.txt file for the metagene_plot_5UTR_CDS_3UTR_v2.c script

final file (a tab-delimited file called “in_features.txt”. No header row.):
# Column 1: transcript name (to match the transcript name in the SAM file). 
# Column 2: total length of the transcript in nucleotides.
# Column 3: position of the first nucleotide of the start codon of the main 
#   Open Reading Frame on the transcript. 
# Column 4: position of the last nucleotide of the stop codon of the main 
#   Open Reading Frame on the transcript.
"""


def get_gene_ID(string):
    name_trunk = string.split('_')
    return name_trunk[0]

def load_UTR(in_path):
    UTR = []
    names = []
    with open(in_path) as file:
        for line in file:
            if len(line) >= 2:
                words = line.split()
                feature = [get_gene_ID(words[0]), words[2], words[3]]
                UTR.append(feature)
                names.append(get_gene_ID(words[0]))

    return UTR, names


def get_3_5_UTR(UTR_3, UTR_5, names_5):
    all_UTR=[]
    counter=0
    for row_3 in UTR_3:
        if row_3[0] in names_5:
            for row_5 in UTR_5:
                if row_3[0] == row_5[0]:
                    counter +=1
                    in_both_lists=[row_3[0], row_5[1], row_5[2], row_3[1], row_3[2]]
                    all_UTR.append(in_both_lists)
    return all_UTR


# wanted: name, total length, pos start, pos end
def get_intervals(all_UTR):
    final_features_all = []
    final_features_meaningfull = []
    counter=0
    for row in all_UTR:
        total_length = int(row[4]) - int(row[1])
        # if forward strand
        if total_length > 0:
            start = int(row[2]) - int(row[1])
            end = int(row[3]) - int(row[1])
        # if reverse strand
        else:
            total_length = int(row[2]) - int(row[3])
            start = total_length - (int(row[1]) - int(row[3]))
            end = total_length - (int(row[4]) - int(row[3]))
        feature = [row[0], total_length, start, end]
        
        final_features_all.append(feature)
    return final_features_all



def print_outputfile(final_features, out_path):
    with open(out_path, 'w', newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')  
        for row in final_features:
            tsv_writer.writerow(row)

def get_all_features(in_path_3, in_path_5, out_final_UTR):
    UTR_3, names_3 = load_UTR(in_path_3)
    UTR_5, names_5 = load_UTR(in_path_5)
    all_UTR = get_3_5_UTR(UTR_3, UTR_5, names_5)
    final_features = get_intervals(all_UTR)
    print_outputfile(final_features, out_final_UTR)
    
    
in_path_3 = "C:/Users/user/Desktop/BPC_FP/programming/5_UTR/3_UTR.gff3"
in_path_5 = "C:/Users/user/Desktop/BPC_FP/programming/5_UTR/5_UTR.gff3"
out_final_UTR = "C:/Users/user/Desktop/BPC_FP/programming/5_UTR/features.txt"
    
get_all_features(in_path_3, in_path_5, out_final_UTR)


        
        











