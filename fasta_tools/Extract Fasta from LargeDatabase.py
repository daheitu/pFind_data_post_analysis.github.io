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
wkpath = r"D:\FastaDatabase\DATABASE_EXTRACT\transcolon"
os.chdir(wkpath)
FastaFile = "uniprot-proteome_UP000006906_con.fasta"
new_fasta = "trans_small.fasta"


def Get_ProList_from_pBuildFile(filename):
    protein_list = []
    pro_simplify_list = []
    f = open(filename, 'r').readlines()
    for line in f[1:]:
        line_list = line.rstrip("\n").split("\t")
        protein_list.append(line_list[0])
        # group = line_list[-2]
        # if group == "":
        #     continue
        # else:
        #     Subset_Pro_List = group.split("/")[:-1]
        #     protein_list += Subset_Pro_List

    pro_simplify_list = []
    for protein in protein_list:
        if protein not in pro_simplify_list:
            pro_simplify_list.append(protein)
    return pro_simplify_list


def Extract_protein(FastaName, tgt_list):
    fasta = open(FastaName, 'r').readlines()
    b = open(new_fasta, "w")
    pro_up_pos = []
    pro_down_pos = []
    protein_total_list = []
    rev_pro = []
    for i in range(len(fasta)):
        if fasta[i][0] == '>':
            pro_up_pos.append(i)
            protein_total_list.append(fasta[i][1:fasta[i].find(" ", 1)])

    for j in range(1, len(pro_up_pos)):
        pro_down_pos.append(pro_up_pos[j] - 1)

    pro_down_pos.append(len(fasta) - 1)
    print(len(pro_up_pos), len(pro_down_pos))

    for pro in tgt_list:
        for k in range(len(protein_total_list)):
            if protein_total_list[k] == pro:
                for i in range(pro_up_pos[k], pro_down_pos[k] + 1):
                    b.write(fasta[i])
                break
            elif k == len(protein_total_list) - 1:
                print(pro)
            else:
                continue
    b.close()
    return


def main():
    # os.chdir(r"E:\workspace\pFindTask86\result\pBuild_tmp")
    FileName_list = os.listdir(os.getcwd())
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
