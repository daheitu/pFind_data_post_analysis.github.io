# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 09:48:50 2016
How to use this script. Copy the pBuild protein report "txt" file, 
Large Database fasta file and this script to one decretory.
cmd 
python Extract Fasta from LargeDatabase.py
@author: Yong Cao
"""

import os
from unicodedata import name
from xml.dom import NamespaceErr

wk_dir = r"Z:\pFind_work_space\pFindTask3\result\pBuild_tmp"
os.chdir(wk_dir)

def get_fasta():
    for fl in os.listdir(wk_dir):
        if fl.endswith(".fasta"):
            return fl
    return "None"

FastaFile = get_fasta()
new_fasta = FastaFile[:-6] + "_filter.fasta"
psm_cuoff = 2


def Get_ProList_from_pBuildFile(filename):
    protein_set = set()
    f = open(filename, 'r').readlines()
    for line in f[1:]:
        line_list = line.rstrip("\n").split("\t")
        protein_name = line_list[1]
        if not protein_name.startswith("REV_"):
            protein_set.add(protein_name)
        
        group = line_list[-2]
        if group == "":
            continue
        else:
            Subset_Pro_List = group.split("/")[:-1]
            for pro in Subset_Pro_List:
                if not pro.startswith("REV_"):
                    protein_set.add(pro)
                
    return protein_set


def Extract_protein(FastaName, name_set):
    fasta = open(FastaName, 'r').readlines()
    b = open(new_fasta, "w")
    i = 0
    writed_pros = set()
    while i < len(fasta) and fasta[i][0] == '>':
        if fasta[i][1:5] == "CON_":
            for m in range(len(fasta[i])):           
                if fasta[i][m] == "|":  # in [" ", "|", "\n", "\t"]:
                    line_Pro_name = fasta[i][1:m]
                    break
                else:
                    continue
        else:
            for m in range(len(fasta[i])):           
                if fasta[i][m] in [" ", "\n", "\t"]:
                    line_Pro_name = fasta[i][1:m]
                    break
                else:
                    continue
        p = i + 1
        if line_Pro_name in name_set:
            # print(line_Pro_name)
            writed_pros.add(line_Pro_name)
            b.write(fasta[i])
            while p < len(fasta) and fasta[p][0] != '>':
                b.write(fasta[p])
                p += 1
        else:
            while p < len(fasta) and fasta[p][0] != '>':
                p += 1
        i = p
    
    over = name_set & writed_pros
    if over == name_set:
        print("well done")
    else:
        print(name_set - writed_pros)
        print(writed_pros - name_set)
    return


def main():
    FileName_list = os.listdir(os.getcwd())
    Total_Pro = []
    for fl in FileName_list:
        if fl.endswith("selected_protein_result.txt"):
            Total_Pro = Get_ProList_from_pBuildFile(fl)
        else:
            continue
    print(len(Total_Pro))
    Extract_protein(FastaFile, Total_Pro)


if __name__ == "__main__":
    main()
