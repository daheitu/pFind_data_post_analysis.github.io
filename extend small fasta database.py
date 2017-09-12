# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 09:48:50 2016

@author: l
"""

import os
os.chdir(r"D:\DATABASE_EXTRACT")
f1=open("uniprot-proteome_hamster_with_virus.fasta")
fastalist = f1.readlines()
namenum = []
for i in range(len(fastalist)):
    if fastalist[i][0]=='>':
        namenum.append(i)
namenum.append(len(fastalist))
print(namenum)

f2=open("list.txt")
namelist=f2.readlines()
f3=open("uniprot-proteome_hamster_with_virus_sub_database.fasta", 'w')
for m in range(len(namelist)):    
    for k in range(len(fastalist)):
        if fastalist[k][1:fastalist[k].find(" ",3)].strip() == namelist[m].strip():
            print(k)
            f3.write(str(fastalist[k].strip()))
            f3.write('\n')
            for n in range(k+1, namenum[namenum.index(k)+1]):
                f3.write(str(fastalist[n].strip()))
            f3.write('\n')
f3.close

        

