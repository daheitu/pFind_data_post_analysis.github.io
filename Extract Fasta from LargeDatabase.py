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
FastaFile = "UniProt-Trembl_human_na_2013-12-5_yhding_con.fasta"
new_fasta = "extract.fasta"


def Get_ProList_from_pBuildFile(filename):
    protein_list = []
    Pro_list_simp = []
    f = open(filename, 'r').readlines()
    for line in f[1:]:
        line_list = line.rstrip("\n").split("\t")
        protein_list.append(line_list[1])
        group = line_list[-2]
        if group == "":
            continue
        else:
            Subset_Pro_List = group.split("/")[:-1]
            protein_list += Subset_Pro_List

    Pro_list_simp = []
    for protein in protein_list:
        if protein not in Pro_list_simp:
            Pro_list_simp.append(protein)
    return Pro_list_simp


def Extract_protein(FastaName, namelist):
    fasta = open(FastaName, 'r').readlines()
    b = open(new_fasta, "w")
    pro_up_pos = []
    pro_down_pos = []
    pro_name_list = []
    rev_pro = []
    for i in range(len(fasta)):
        if fasta[i][0] == '>':
            pro_up_pos.append(i)
            pro_name_list.append(fasta[i][1:fasta[i].find(" ", 1)])

    for j in range(1, len(pro_up_pos)):
        pro_down_pos.append(pro_up_pos[j] - 1)

    pro_down_pos.append(len(fasta) - 1)
    print(len(pro_up_pos), len(pro_down_pos))

    for pro in namelist:
        for k in range(len(pro_name_list)):
            if pro_name_list[k] == pro:
                for i in range(pro_up_pos[k], pro_down_pos[k] + 1):
                    b.write(fasta[i])
                break
            elif k == len(pro_name_list) - 1:
                print(pro)
            else:
                continue
    b.close()
    return


def main():
    os.chdir(r"C:\Users\Yong\Desktop\sgc")
    FileName_list = os.listdir(os.getcwd())
    fastfile = ""
    Total_Pro = []
    Total_Pro_simp = []
    for fl in FileName_list:
        if fl[-4:] == ".txt":
            Total_Pro += Get_ProList_from_pBuildFile(fl)
        else:
            continue

    for pro in Total_Pro:
        if pro not in Total_Pro_simp:
            Total_Pro_simp.append(pro)
    print(len(Total_Pro_simp))
    Extract_protein(FastaFile, Total_Pro_simp)


if __name__ == "__main__":
    main()
