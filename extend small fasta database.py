# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 09:48:50 2016

@author: Yong Cao
"""

import os
os.chdir(r"D:\fasta database\DATABASE_EXTRACT")  # work directory
f1 = open("UniProt-Trembl_human_na_2013-12-5_yhding_con.fasta"
          )  # your orignal fasta file to extrac small
fastalist = f1.readlines()
pro_up_pos = []
pro_down_pos =[]
pro_name_list = []
for i in range(len(fastalist)):
    if fastalist[i][0] == '>':
        pro_up_pos.append(i)
        pro_name_list.append(fastalist[i][1:fastalist[i].find(" ", 1)])
        
for j in range(1, len(pro_up_pos)):
    pro_down_pos.append(pro_up_pos[j]-1)

pro_down_pos.append(len(fastalist)-1)
print(len(pro_up_pos), len(pro_down_pos))

f2 = open("list.txt")  # protein name list for your requir
namelist = f2.readlines()
f3 = open("Jingjidatabase.fasta", 'w')
for m in range(len(namelist)):
    for k in range(len(pro_name_list)):
        if pro_name_list[k] == namelist[m].rstrip("\n"):
            for i in range(pro_up_pos[k], pro_down_pos[k]+1):
                f3.write(fastalist[i])
            break
        elif k == len(pro_name_list)-1 :
            print(namelist[m])
f3.close()
