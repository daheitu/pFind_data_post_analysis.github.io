# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:38:11 2016
@author: CaoYong
"""

import os
os.chdir(r"D:\softwareData\plink2\search_task_2017.09.05.15.20.16_BSA_low_concentration\reports")
filename = os.listdir(os.getcwd())
for fl in filename:
    if fl[-22:-4] == "cross-linked_sites":
        f=open(fl)
        site_table=f.readlines()
        ordr_pos=[]
        sit_list=[]
        Peptide=[]
        for i in range(2,len(site_table)):
            if len([x for x in site_table[i].strip().split(',') if x !=''])==4:
                ordr_pos.append(i)
        ordr_pos.append(len(site_table))
        colum_name=[]
        for m in range(len(ordr_pos)-1):
            for k in range(ordr_pos[m]+1,ordr_pos[m+1]):
                a=site_table[k].strip().split(',')[2].find(".")
                colum_name.append(site_table[k].strip().split(',')[2][:a])
        colum_name=list(set(colum_name))
        #print colum_name
        col=["linked_site","total_spec","best_svm_score","peptide","peptide_mass","prote_type"]
        for i in range(len(colum_name)):
            col.append(str(colum_name[i])+"_spec")
            col.append(str(colum_name[i])+"_svm")
            col.append(str(colum_name[i])+"unique_pep_num")
        b=open("report.txt", 'w')
        b.write('\t'.join(col))
        b.write('\n')
        m=0
        for l in range(len(ordr_pos)-1):
            vars()["sub_linked_if"+str(l)]=[]
            vars()["sub_linked_if"+str(l)].append(site_table[ordr_pos[l]].strip().split(',')[1])
            vars()["sub_linked_if"+str(l)].append(site_table[ordr_pos[l]].strip().split(',')[3])
            vars()["sub_linked_if"+str(l)].append(site_table[ordr_pos[l]+1].strip().split(',')[9])
            vars()["sub_linked_if"+str(l)].append(site_table[ordr_pos[l]+1].strip().split(',')[5])
            vars()["sub_linked_if"+str(l)].append(site_table[ordr_pos[l]+1].strip().split(',')[4])
            linked_site=site_table[ordr_pos[l]].strip().split(',')[1]
            m=linked_site.find("(")
            n=linked_site.find(")")
            p=linked_site.find("-")
            x=linked_site.find("(",p)
            y=linked_site.find(")",p)
            if linked_site[:m]==linked_site[p+1:x]:
                vars()["sub_linked_if"+str(l)].append("intra")
            else:
                vars()["sub_linked_if"+str(l)].append("inter")
            colum_spec=[]
            colum_svm=[]
            raw_pep=[]
            for h in range(len(colum_name)):
                vars()[colum_name[h]+"spec"]=0
                vars()[colum_name[h]+"svm"]=[]
                vars()[colum_name[h]+"pep"]=[]
                for j in range(ordr_pos[l]+1,ordr_pos[l+1]):            
                    if colum_name.index(site_table[j].strip().split(',')[2][:site_table[j].strip().split(',')[2].find(".")])==h:
                        vars()[colum_name[h]+"spec"]+=1
                        vars()[colum_name[h]+"svm"].append(site_table[j].strip().split(',')[9])
                        vars()[colum_name[h]+"pep"].append(site_table[j].strip().split(',')[5])
                if vars()[colum_name[h]+"spec"]==0:
                    vars()["sub_linked_if"+str(l)].append("")
                    vars()["sub_linked_if"+str(l)].append("")
                    vars()["sub_linked_if"+str(l)].append("")
                else:
                    vars()["sub_linked_if"+str(l)].append(str(vars()[colum_name[h]+"spec"]))
                    vars()["sub_linked_if"+str(l)].append(vars()[colum_name[h]+"svm"][0])
                    vars()["sub_linked_if"+str(l)].append(str(len(list(set(vars()[colum_name[h]+"pep"])))))
            #print vars()["sub_linked_if"+str(l)]
            
            b.write('\t'.join(vars()["sub_linked_if"+str(l)]))
            b.write('\n')
        b.close()

"""
            b = open("report.txt", 'a+')
            b.seek(0)
            rep_table = b.readlines()
            col_dic = {}
            total_spectra = 0
            total_colom = len(rep_table[0].strip().split("\t"))
            intra_num = 0
            distance_list = []
            for i in range(1, len(rep_table)):
                total_spectra += int(rep_table[i].strip("\n").split("\t")[1])
                if rep_table[i].strip("\n").split("\t")[5] == "intra":
                    intra_num += 1
                if rep_table[i].strip("\n").split("\t")[6]:
                    distance_list.append(rep_table[i].strip("\n").split("\t")[6])

            col_dic[5] = float(intra_num) / (len(rep_table) - 1)
            col_dic[6] = ""
            col_dic[0] = len(rep_table) - 1
            col_dic[1] = total_spectra
            col_dic[2] = ""
            col_dic[3] = ""
            col_dic[4] = ""

            for k in [7, 9]:
                for j in range(k, total_colom, 3):
                    vars()["spec" + str(j) + "_list"] = []
                    for i in range(1, len(rep_table)):
                        # print rep_table[i].strip("\n").split("\t")
                        if rep_table[i].strip("\n").split("\t")[j]:
                            vars()["spec" + str(j) + "_list"].append(int(rep_table[i].strip("\n").split("\t")[j]))
                    col_dic[j] = sum(vars()["spec" + str(j) + "_list"])

            for j in range(8, total_colom, 3):
                vars()["col_svm" + str(j) + "_list"] = []
                for i in range(1, len(rep_table)):
                    if rep_table[i].strip("\n").split("\t")[j]:
                        vars()["col_svm" + str(j) + "_list"].append(
                            round(float(rep_table[i].strip("\n").split("\t")[j]), 2))
                col_dic[j] = str(min(vars()["col_svm" + str(j) + "_list"])) + ',' + str(
                    max(vars()["col_svm" + str(j) + "_list"]))

            last = []
            for k in range(total_colom):
                last.append(str(col_dic[k]))
            # print col_dic
            b.write("\t".join(last))
            b.close()

"""
