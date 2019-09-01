# coding: utf-8
"""
This script can help you to summary the plink2 report file
"""

import os
import re 

os.chdir(
    r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.08.26.13.00.31\reports"
)
spec_cutoff = 0  # spectra number cut-off
Best_evalue_cutoff = 2
E_value_cutoff_SpecLvl = 2


def judgeHomoHetro(linked_site, pepXL):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]

    if protein1 != protein2:
        return "Inter"
    else:
        linkPosPep1, linkPosPep2 = re.findall("\((\d*)\)", pepXL)
        p1 = pepXL.find(")-")
        m1 = pepXL.find("(" + linkPosPep1 + ")")
        n1 = pepXL.find("(" + linkPosPep2 + ")", p1)
        pep1 = pepXL[:m1]
        pep2 = pepXL[p1 + 2: n1]
        deltaList = [int(position1) - int(linkPosPep1), int(position2) - int(linkPosPep2)]
        pep1IDX = [1+deltaList[0], len(pep1)+ deltaList[0]]
        pep2IDX = [1+deltaList[1], len(pep2)+ deltaList[1]]
        if pep1IDX[1] < pep2IDX[0] or pep1IDX[0] > pep2IDX[1]:
            return "Intra"
        else:
            return "Inter"


def site_list_process(site_list, pepXLlist):
    i = 0
    while i < len(site_list):
        if "REVERSE" in site_list[i] or "gi|CON" in site_list[i]:
            site_list.remove(site_list[i])
        else:
            i += 1

    if site_list == "":
        return "", None
    else:
        link_type_list = []
        for i in range(len(site_list)):
            #boolInter = False
            for j in range(len(pepXLlist)):
                if judgeHomoHetro(site_list[i], pepXLlist[j]) == "Inter":
                    link_type_list.append("Inter")
                    break
                else:
                    continue
            if j == len(pepXLlist)-1:
                link_type_list.append("Intra")
        linkType = "/".join(link_type_list)
        site = "/".join(site_list)
        return site, linkType
    

def get_report_file_name():
    path = os.getcwd()
    path_d = os.path.dirname(path)
    os.chdir(path_d)
    
    file_list = os.listdir(path_d)
    for fl in file_list:
        if fl[-5:] in ["pfind", "plink"]:
            para = open(fl).readlines()
            for line in para:
                #print(line)
                if line[:10] == "spec_title":
                    spec_title = line.rstrip("\n").split("=")[1].strip()
                if line[:11] == "enzyme_name":
                    enzyme = line.rstrip("\n").split("=")[1].strip()
                if line[:7] == "linker1":
                    linker = line.rstrip("\n").split("=")[1].strip()
        else:
            continue

    report_file_name = spec_title + "_" + enzyme + "_" + linker + "v5.txt"
    os.chdir(path)
    return report_file_name


def get_crosslink_site_info(site_table, b):
    raw_name_list = []
    for line in site_table[2:]:
        line_list = line.rstrip("\n").split(",")
        if line_list[0] == "":
            raw_name = line_list[2][:line_list[2].find(".")]
            if raw_name not in raw_name_list:
                raw_name_list.append(raw_name)
    raw_name_list.sort()
    print("All raw files are: " + ",".join(raw_name_list))
    col = [
        "Linked Site", "Total Spec", "Best E-value", "Best Svm Score",
        "Peptide", "Inter or Intra Molecular"
    ]
    for name in raw_name_list:
        col.append(name + "_SpecNum")
        col.append(name + "_E-value")
        col.append(name + "_UniquePepNum")

    b.write('\t'.join(col))
    b.write('\n')
    return raw_name_list


def cal_sumOfOneColumn(openedfl, k_column):
    f = openedfl; k = k_column
    sumNum = 0
    for line in f[1:]:
        lineList = line.rstrip("\n").split("\t")
        val = lineList[k]
        if val:
            sumNum += int(lineList[k])
    return sumNum


def cal_numRange(openedfl, k_column):
    f = openedfl; k = k_column
    valList = []
    for line in f[1:]:
        lineList = line.rstrip("\n").split("\t")
        val = lineList[k]
        if val:
            valList.append(val)
    newList = sorted(valList, key = lambda x: float(x))
    return newList[0] + "--" + newList[-1]
    

def statistic_report_file():
    report_file_name = get_report_file_name()
    c = open(report_file_name, 'a')
    rep_table = open(report_file_name, 'r').readlines()
    title_list = rep_table[0].split("\t")
    col_dic = {}
    total_colom = len(rep_table[0].strip().split("\t"))
    intra_num = 0
    for i in range(1, len(rep_table)):
        if rep_table[i].strip("\n").split("\t")[5] == "Intra":
            intra_num += 1
    col_dic[5] = float(intra_num) / (float(len(rep_table)) - 1.0)
    col_dic[0] = len(rep_table) - 1
    col_dic[1] = cal_sumOfOneColumn(rep_table, 1)
    col_dic[2] = ""
    col_dic[3] = ""
    col_dic[4] = ""

    column_sub_dic = {}
    for k in [6, 8]:
        for j in range(k, total_colom, 3):
            #print("the column is " + str(j))
            col_dic[j] = cal_sumOfOneColumn(rep_table, j)

    for j in range(7, total_colom, 3):
        col_dic[j] = cal_numRange(rep_table, j)
        
    last = []
    for k in range(total_colom):
        last.append(str(col_dic[k]))
    c.write("\t".join([str(ele) for ele in last]) + "\n")
    raw_dic = {}
    for i in range(7, len(title_list)):
        raw_name_list = title_list[i].split("_")[:-1]
        raw_name = "_".join(raw_name_list)
        if raw_name not in raw_dic:
            raw_dic[raw_name] = [raw_name, last[i]]
        else:
            raw_dic[raw_name].append(last[i])

    raw_list = list(raw_dic.keys())
    raw_list.sort()
    for raw in raw_list:
        c.write("\t".join(raw_dic[raw]))
        c.write("\n")

    c.close()
    return


def splitResult(openedfl, raw_name_list):
    finalList = []
    f = openedfl
    n = 2
    while n < len(f):
        spec_dic = {}
        pep_dic = {}
        evalue_dic = {}

        line = f[n]
        line_list = line.rstrip("\n").split(",")
        if line_list[0].isdigit():
            site_list = [line_list[1]]
            p = n + 1
        else:
            print(n)
        
        cell1 = f[p].rstrip("\n").split(",")[0]
        
        if cell1.isdigit():
            pass
        else:
            while p < len(f):
                line_list = f[p].rstrip("\n").split(",")
                if line_list[0] == "" and line_list[1].isdigit():
                    break
                else:
                    if line_list[0] == "SameSet":
                        site_list.append(line_list[1])              
                    p += 1
        
        site = site_list_process(site_list, ["WFC(2)-XSV(2)"])[0]
        if site == "":
            while p < len(f):
                if f[p].rstrip("\n").split(",")[0].isdigit():
                    break
                else:
                    p += 1
        else:
            bestSVMscore = f[p].rstrip("\n").split(",")[9]
            pep_std_list = [] # f[p].rstrip("\n").split(",")[5]
            while p < len(f): # and f[p].rstrip("\n").split(",")[0] == "":
                line_list = f[p].rstrip("\n").split(",")
                if line_list[0].isdigit():
                    break
                else:
                    pep_std = line_list[5]
                    if pep_std not in pep_std_list:
                        pep_std_list.append(pep_std)
        
                    evalue = line_list[8]
                    if float(evalue) < E_value_cutoff_SpecLvl:
                        pep = line_list[5]
                        raw_name = line_list[2][:line_list[2].find(".")]                        
                        if raw_name not in spec_dic:
                            spec_dic[raw_name] = 1
                        else:
                            spec_dic[raw_name] += 1

                        if raw_name not in pep_dic: 
                            pep_dic[raw_name] = [pep]
                        else:
                            if pep not in pep_dic[raw_name]:
                                pep_dic[raw_name].append(pep)

                        if raw_name not in evalue_dic:
                            evalue_dic[raw_name] = evalue
                        else:
                            if float(evalue) < float(evalue_dic[raw_name]):
                                evalue_dic[raw_name] = evalue

                    p += 1

            totalPep = ";".join(pep_std_list)
            link_type = site_list_process(site_list, pep_std_list)[1]
            total_spec = 0
            min_evalue = 1
            for key in evalue_dic:
                total_spec += spec_dic[key]
                if float(evalue_dic[key]) < float(min_evalue):
                    min_evalue = evalue_dic[key]
                else:
                    continue
            if total_spec > spec_cutoff and float(min_evalue) < Best_evalue_cutoff:
                rep_list = [site, total_spec, min_evalue,\
                    bestSVMscore, totalPep, link_type]
                
                for raw in raw_name_list:
                    if raw not in spec_dic:
                        SEP = ["", "", ""]
                    else:
                        SEP = [spec_dic[raw], evalue_dic[raw], len(pep_dic[raw])]
                    rep_list.extend(SEP)

                finalList.append("\t".join([str(ele) for ele in rep_list]))

        n = p
    return finalList


def main():
    filename_list = os.listdir(os.getcwd())
    link_site_file = ""
    for fl in filename_list:
        if fl[-22:-4] == "cross-linked_sites":
            link_site_file = fl
        else:
            continue
    if link_site_file == "":
        print("Please check your cross-linked_sites file")
    else:
        f = open(link_site_file, 'r').readlines()

    report_file_name = get_report_file_name()
    b = open(report_file_name, 'w')
    raw_name_list = get_crosslink_site_info(f, b)
    finalList = splitResult(f, raw_name_list)
    finalList = sorted(finalList, key = lambda x: int(x.split("\t")[1]), reverse = True)
    for line in finalList:
        b.write(line + "\n")
    b.close()
    print("Well Done")
    #statistic_report_file()


if __name__ == "__main__":
    main()