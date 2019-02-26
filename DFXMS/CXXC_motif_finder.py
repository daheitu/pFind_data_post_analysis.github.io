# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 21:47:36 2018

@author: Yong
"""

import os
import re
from bs4 import BeautifulSoup
import urllib3
import time
import random

database_path = r"D:\fasta database\uniprot_Ecoli_K12.fasta"

mtf = "C[GAVLISTKRHDNEQPWFYM]{6,}C[GAVLISTHDNEQPWFYM]{0,6}C[GAVLISTKRHDNEQPWFYM]{6,}C"
mf2 = "C[GAVLISTHDNEQPWFYM]{0,6}C"
# re.split()

"""
seq = ">BSA \n MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQ\
CPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"
"""


def motif_finder(seq):
    mtf_list = re.findall(mtf, seq)
    return mtf_list




def get_SubLocation(Uniprot_ID):
    url = "https://www.uniprot.org/uniprot/" + Uniprot_ID + ".txt"
    urllib3.disable_warnings()
    # url = "https://www.uniprot.org/uniprot/P0A953.txt"

    http = urllib3.PoolManager()
    respon = http.request('GET', url)
    soup = BeautifulSoup(respon.data, "html5lib")
    soup_txt = str(soup).split("\n")
    sub_loc = ""
    for line in soup_txt:
        # print(line)
        if line[:29] == "CC   -!- SUBCELLULAR LOCATION":
            sub_loc = line[30:-1] 
            print(sub_loc)
            break
        else:
            continue
    return sub_loc


fasta = open(database_path, 'r').readlines()
b = open("report.txt", 'w')
sleep_time = random.uniform(1.2, 5)
i = 0
while i < len(fasta) and fasta[i][0] == '>':
    # print(i)
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
    
    pro_seq = ""
    while p < len(fasta):
        if fasta[p][0] == '>':
            break
        else:
            pro_seq += fasta[p].strip()
        p += 1
    mtf_list = motif_finder(pro_seq)
    mtf_flt_list = []
    for word in mtf_list:
        motif = re.search(mf2, word).group()
        part_list = word.split(motif)
        dgst_bool = True
        for part in part_list:
            if "K" not in part and "R" not in part:
                dgst_bool = False
                break
            else:
                continue
        if dgst_bool == True:
            mtf_flt_list.append(word)
        else:
            continue
    w_list = [line_Pro_name]
    if mtf_flt_list:
        uni_id = line_Pro_name.split("|")[1]
        location = get_SubLocation(uni_id)
        w_list.extend(mtf_flt_list)
        w_list.append(location)
        b.write("\t".join(w_list) + "\n")
    else:
        pass
    print(p)
    i = p 
    time.sleep(sleep_time)
b.close()
