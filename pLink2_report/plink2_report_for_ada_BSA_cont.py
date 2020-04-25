# coding: utf-8
"""
This script can help you to summary the plink2 report file
"""

import os
import re
from calSynthiticNum import calsyntheticNum, calsyntheticNumforFix

os.chdir(
    r"G:\msData\synthetic_pepteide_rawdata\CV1_CV3\DSS\CV1_CV3_DSS\reports"
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
        pep2 = pepXL[p1 + 2:n1]
        deltaList = [
            int(position1) - int(linkPosPep1),
            int(position2) - int(linkPosPep2)
        ]
        pep1IDX = [1 + deltaList[0], len(pep1) + deltaList[0]]
        pep2IDX = [1 + deltaList[1], len(pep2) + deltaList[1]]
        if pep1IDX[1] < pep2IDX[0] or pep1IDX[0] > pep2IDX[1]:
            return "Intra"
        else:
            return "Inter"


def site_list_process(site_list, pepXLlist):
    con_site_list = []
    i = 0
    while i < len(site_list):
        if "gi|CON" in site_list[i] or "sp|" in site_list[i]:
            list_pop = site_list.pop(i)
            #print(list_pop)
            con_site_list.append(list_pop)
        else:
            i += 1
    print(con_site_list)
    print(site_list)
    if site_list == []:
        return "contam", "/".join(con_site_list), "Inter"
    else:
        link_type_list = []
        for i in range(len(site_list)):
            boolInter = False
            for j in range(len(pepXLlist)):
                if judgeHomoHetro(site_list[i], pepXLlist[j]) == "Inter":
                    boolInter = True
                    break
                else:
                    continue
            if boolInter:
                link_type_list.append("Inter")
            else:
                link_type_list.append("Intra")
        linkType = "/".join(link_type_list)
        site = "/".join(site_list)
        return "BSA", site, linkType


def get_report_file_name(): 
    path_d = os.path.dirname(os.getcwd())
    file_list = os.listdir(path_d)
    for fl in file_list:
        if fl.endswith("plink"):
            para = open(os.path.join(path_d, fl)).readlines()
            for line in para:
                if line.startswith("spec_title"):
                    spec_title = line.split("=")[1].strip()
                if line.startswith("enzyme_name"):
                    enzyme = line.split("=")[1].strip()
                if line.startswith("linker1"):
                    linker = line.split("=")[1].strip()
        else:
            continue

    report_file_name = spec_title + "_" + enzyme + "_" + linker + "_v5.csv"
    return report_file_name


def get_crosslink_site_info(site_table):
    raw_name_list = []
    for line in site_table[2:]:
        line_list = line.rstrip("\n").split(",")
        if line_list[0] == "":
            raw_name = line_list[2][:line_list[2].find(".")]
            if raw_name not in raw_name_list:
                raw_name_list.append(raw_name)
    raw_name_list.sort()
    print("All raw files are: " + ",".join(raw_name_list))
    

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
            #print("the column is " + str(j))
            col_dic[j] = cal_sumOfOneColumn(rep_table, j)

    for j in range(7, total_colom, 3):
        col_dic[j] = cal_numRange(rep_table, j)

    last = []
    for k in range(total_colom):
        last.append(str(col_dic[k]))
    c.write(",".join([str(ele) for ele in last]) + "\n")
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
        c.write(",".join(raw_dic[raw]))
        c.write("\n")

    c.close()
    return


def add_info2repList(openedfl, startPos, spec_dic, pep_dic, evalue_dic, link_type, site, raw_name_list, finalList):
    f = openedfl
    p = startPos
    bestSVMscore = f[p].rstrip("\n").split(",")[9]
    pep_std_list = []  
    while p < len(f):  
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
        rep_list = [site, total_spec, min_evalue, bestSVMscore, totalPep, link_type]

        for raw in raw_name_list:
            if raw not in spec_dic:
                SEP = ["", "", ""]
            else:
                SEP = [
                    spec_dic[raw], evalue_dic[raw],
                    len(pep_dic[raw])
                ]
            rep_list.extend(SEP)

        finalList.append(",".join([str(ele) for ele in rep_list]))
    return p



def splitResult(openedfl, raw_name_list):
    con_finalList = []
    bsa_finalList = []
    f = openedfl
    n = 2
    while n < len(f):

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
        print(site_list)
        tp_bool, site, link_type = site_list_process(site_list, ["WFC(2)-XSV(2)"])
        print(tp_bool, site, link_type)
        spec_dic = {}
        pep_dic = {}
        evalue_dic = {}
        if tp_bool == "contam":
            end_pos = add_info2repList(f, p, spec_dic, pep_dic, evalue_dic, link_type, site, raw_name_list, con_finalList)
        else:
            end_pos = add_info2repList(f, p, spec_dic, pep_dic, evalue_dic, link_type, site, raw_name_list, bsa_finalList)

        n = end_pos
    return con_finalList, bsa_finalList


def wirte_repList2FL(repList, flname, raw_name_list):
    repList = sorted(repList, key=lambda x: int(x.split(",")[1]),
                       reverse=True)
    b = open(flname, 'w')
    col = [
        "Linked Site", "Total Spec", "Best E-value", "Best Svm Score",
        "Peptide", "Inter or Intra Molecular"
    ]
    for name in raw_name_list:
        col.append(name + "_SpecNum")
        col.append(name + "_E-value")
        col.append(name + "_UniquePepNum")
    b.write(",".join(col)+"\n")
    for line in repList:
        b.write(line + "\n")
    b.close()


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

    raw_name_list = get_crosslink_site_info(f)
    con_finalList, bsa_finalList = splitResult(f, raw_name_list)
    wirte_repList2FL(bsa_finalList, report_file_name, raw_name_list)
    wirte_repList2FL(con_finalList, "con_reprter.csv", raw_name_list)
    print("Well Done")
    statistic_report_file(report_file_name)
    statistic_report_file("con_reprter.csv")


if __name__ == "__main__":
    main()