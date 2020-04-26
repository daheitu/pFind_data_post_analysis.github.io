# coding: utf-8
"""
This script can help you to summary the plink2 report file
"""

import os
import re
from calSynthiticNum import calsyntheticNum, calsyntheticNumforFix

os.chdir(
    r"F:\Script\test\ada"
)
spec_cutoff =   0# 大于spectra number cut-off
Best_evalue_cutoff = 2
E_value_cutoff_SpecLvl = 2


def judge_bsa_type(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]
    if protein1 == "BSA" and protein2 == "BSA":
        return True
    else:
        return False


def site_list_process(site_list):
    con_site_list = []
    i = 0
    while i < len(site_list):
        if judge_bsa_type(site_list[i]):
            i += 1
        else:
            con_site_list.append(site_list.pop(i))
    if site_list:
        return True, site_list
    else:
        return False, con_site_list


def get_crosslink_site_info(site_table):
    raw_name_list = []
    for line in site_table[2:]:
        line_list = line.rstrip("\n").split(",")
        if line_list[0] == "":
            raw_name = line_list[2][:line_list[2].find(".")]
            if raw_name not in raw_name_list:
                raw_name_list.append(raw_name)
    raw_name_list.sort()
    return raw_name_list


def cal_sumOfOneColumn(openedfl, k_column):
    f = openedfl
    k = k_column
    sumNum = 0
    for line in f[1:]:
        lineList = line.rstrip("\n").split(",")
        val = lineList[k]
        if val:
            sumNum += int(lineList[k])
    return sumNum


def cal_numRange(openedfl, k_column):
    f = openedfl
    k = k_column
    valList = []
    for line in f[1:]:
        lineList = line.rstrip("\n").split(",")
        val = lineList[k]
        if val:
            valList.append(val)
    newList = sorted(valList, key=lambda x: float(x))
    return newList[0] + "--" + newList[-1]


def info_summary(site, openedfl, startPos, raw_name_list):
    f = openedfl
    p = startPos
    spec_dic = {}
    pep_dic = {}
    evalue_dic = {}
    bestSVMscore = f[p].rstrip("\n").split(",")[9]
    pep_std_list = []  # f[p].rstrip("\n").split(",")[5]
    while p < len(f):  # and f[p].rstrip("\n").split(",")[0] == "":
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
    link_type = "intra"
    total_spec = 0
    min_evalue = 1
    for key in evalue_dic:
        total_spec += spec_dic[key]
        if float(evalue_dic[key]) < float(min_evalue):
            min_evalue = evalue_dic[key]
        else:
            continue
    if total_spec > spec_cutoff and float(
            min_evalue) < Best_evalue_cutoff:
        rep_list = [site, total_spec, min_evalue,\
            bestSVMscore, totalPep, link_type]

        for raw in raw_name_list:
            if raw not in spec_dic:
                SEP = ["", "", ""]
            else:
                SEP = [
                    spec_dic[raw], evalue_dic[raw],
                    len(pep_dic[raw])
                ]
            rep_list.extend(SEP)

        return rep_list, p
    else:
        return False, p



def splitResult(openedfl, raw_name_list):
    finalList_BSA = []
    finalList_Con = []
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

        while p < len(f):
            line_list = f[p].rstrip("\n").split(",")
            if line_list[0] == "" and line_list[1].isdigit():
                break
            else:
                if line_list[0] == "SameSet":
                    site_list.append(line_list[1])
                p += 1

        bool_isBSA, site_list = site_list_process(site_list)
        site = "/".join(site_list)
        rep_list, n = info_summary(site, f, p, raw_name_list)
        if rep_list and bool_isBSA:
            finalList_BSA.append(",".join([str(ele) for ele in rep_list]))
        elif rep_list and bool_isBSA == False:
            finalList_Con.append(",".join([str(ele) for ele in rep_list]))      
    return finalList_BSA, finalList_Con


def writeFinal2file(finalList, flName, raw_name_list):
    b = open(flName, 'w')
    col = [
        "Linked Site", "Total Spec", "Best E-value", "Best Svm Score",
        "Peptide", "Inter or Intra Molecular"
    ]
    for name in raw_name_list:
        col.append(name + "_SpecNum")
        col.append(name + "_E-value")
        col.append(name + "_UniquePepNum")

    b.write(','.join(col))
    b.write('\n')
    finalList = sorted(finalList,
                       key=lambda x: int(x.split(",")[1]),
                       reverse=True)
    for line in finalList:
        b.write(line + "\n")
    b.close()


def statistic_report_file(report_file_name):
    c = open(report_file_name, 'a')
    rep_table = open(report_file_name, 'r').readlines()
    title_list = rep_table[0].split(",")
    col_dic = {}
    total_colom = len(rep_table[0].strip().split(","))
    intra_num = 0
    for i in range(1, len(rep_table)):
        if rep_table[i].strip("\n").split(",")[5] == "Intra":
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
            col_dic[j] = cal_sumOfOneColumn(rep_table, j)

    for j in range(7, total_colom, 3):
        col_dic[j] = cal_numRange(rep_table, j)

    last = []
    for k in range(total_colom):
        last.append(str(col_dic[k]))
    c.write(",".join([str(ele) for ele in last]) + "\n")
    c.close()



def main():
    filename_list = os.listdir(os.getcwd())
    link_site_file = ""
    for fl in filename_list:
        if fl.endswith("cross-linked_sites.csv"):
            link_site_file = fl
        else:
            continue
    if link_site_file == "":
        print("Please check your cross-linked_sites file")
    else:
        f = open(link_site_file, 'r').readlines()

    raw_name_list = get_crosslink_site_info(f)
    finalList_BSA, finalList_Con = splitResult(f, raw_name_list)
    print(finalList_BSA)
    writeFinal2file(finalList_BSA, "bsa_report.csv", raw_name_list)
    writeFinal2file(finalList_Con, "con_report.csv", raw_name_list)
    statistic_report_file("bsa_report.csv")
    statistic_report_file("con_report.csv")
    print("Well Done")

if __name__ == "__main__":
    main()