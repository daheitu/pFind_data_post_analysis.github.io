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
os.chdir(r"C:\Users\Yong\Desktop\lumos")
FastaFile = "human_unitprot_contaminant.fasta"
new_fasta = "Atg16.fasta"


def Get_ProList_from_pBuildFile(filename):
    protein_list = []
    f = open(filename, 'r').readlines()
    for line in f[1:]:
        line_list = line.rstrip("\n").split("\t")
        protein_name = line_list[1]
        if protein_name not in protein_list:
            protein_list.append(line_list[1])
        else:
            pass
        
        group = line_list[-2]
        if group == "":
            continue
        else:
            Subset_Pro_List = group.split("/")[:-1]
            for pro in Subset_Pro_List:
                if pro not in protein_list:
                    protein_list.append(pro)
                else:
                    pass
    return protein_list


def Extract_protein(FastaName, namelist):
    fasta = open(FastaName, 'r').readlines()
    b = open(new_fasta, "w")
    i = 0
    pro_wr_num = 0
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
        if line_Pro_name in namelist:
            # print(line_Pro_name)
            pro_wr_num += 1
            b.write(fasta[i])
            while p < len(fasta) and fasta[p][0] != '>':
                b.write(fasta[p])
                p += 1
        else:
            while p < len(fasta) and fasta[p][0] != '>':
                p += 1
        i = p
    
    if len(namelist) != pro_wr_num:
        print("missing")
        print(len(namelist), pro_wr_num)
    else:
        print("well done")
    return


def main():
    FileName_list = os.listdir(os.getcwd())
    Total_Pro = []
    for fl in FileName_list:
        if fl[-4:] == ".txt":
            Total_Pro = Get_ProList_from_pBuildFile(fl)
        else:
            continue
    print(len(Total_Pro))
    Extract_protein(FastaFile, Total_Pro)


if __name__ == "__main__":
    main()
