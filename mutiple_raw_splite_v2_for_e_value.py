# coding: UTF-8
"""
Created on Mon Aug 15 10:38:11 2017
@author: CaoYong
This script can help you to deal with the plink2 report file 
"""

import os
import re

os.chdir(
    r"C:\Users\Yong\Documents\pLink\search_task_2018.01.23.11.09.05_SAGA\reports"
)
spec_cutoff = 2  # spectra number cut-off
E_value_cutoff = 0.01


def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find("-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:n].strip()
    return protein1, protein2, position1, position2


def site_correct(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    if "/" in linked_site:
        return linked_site
    else:
        if int(position1) <= int(position2):
            return linked_site
        else:
            a = linked_site.split("-")[0]
            b = linked_site.split("-")[1]
            return b + "-" + a


# 找到包含site以及对应的起始和终止行
def get_site_lines(site_table):
    site_pos_list = []
    site_begain_pos_list = []
    site_end_pos_list = []
    site_pos_dic = {}
    for i in range(2, len(site_table)):
        line_list = site_table[i].rstrip("\n").split(",")
        if line_list[0].isdigit():
            site_pos_list.append(i)
    for pos in site_pos_list:
        site_begain_pos_list.append(pos + 1)
    for pos in site_pos_list[1:]:
        site_end_pos_list.append(pos - 1)

    site_end_pos_list.append(len(site_table) - 1)
    for i in range(len(site_pos_list)):
        site = site_table[site_pos_list[i]].strip().split(",")[1]
        site_beg_end = [site_begain_pos_list[i], site_end_pos_list[i]]
        site_pos_dic[site] = site_beg_end
    return site_pos_dic


# 去掉sameset
def remove_sameset(tab1):
    site_pos_orig_dic = get_site_lines(tab1)
    tmp = open("tmpfile", "w")
    tmp.write(tab1[0])
    tmp.write(tab1[1])
    for site in site_pos_orig_dic:
        [site_up, site_down] = site_pos_orig_dic[site]
        line_up1 = tab1[site_up].rstrip("\n").split(",")[0]
        if line_up1 == "SameSet":
            site_modify = [site]
            for i in range(site_up, site_down + 1):
                line_list = tab1[i].rstrip("\n").split(",")
                first_cell = line_list[0]
                if first_cell == "SameSet":
                    site_modify.append(line_list[1])
                else:
                    continue
            site_line_list = tab1[site_up - 1].rstrip("\n").split(",")
            site_line_list[1] = "/".join(site_modify)
            tmp.write(",".join(site_line_list))
            tmp.write("\n")
            for i in range(site_up, site_down + 1):
                line_list = tab1[i].rstrip("\n").split(",")
                first_cell = line_list[0]
                if first_cell != "SameSet":
                    tmp.write(tab1[i])
                else:
                    continue

        elif line_up1 == "":
            for i in range(site_up - 1, site_down + 1):
                tmp.write(tab1[i])

    tmp.close()
    return


def get_crosslink_site_info(site_table, site_pos_dic):
    b = open("report.txt", 'w')
    raw_name_list = []
    for line in site_table[2:]:
        line_list = line.rstrip("\n").split(",")
        if line_list[0] == "":
            raw_name = line_list[2][:line_list[2].find(".")]
            if raw_name not in raw_name_list:
                raw_name_list.append(raw_name)
    raw_name_list.sort()
    print("All raw file are: " + ",".join(raw_name_list))
    col = [
        "Linked Site", "Total Spec", "Best E-value", "Best Svm Score",
        "Peptide", "Peptide mass", "Prote type"
    ]
    for name in raw_name_list:
        col.append(name + "_SpecNum")
        col.append(name + "_E-value")
        col.append(name + "_UniquePepNum")

    b.write('\t'.join(col))
    b.write('\n')

    for site in site_pos_dic:
        [site_up, site_down] = site_pos_dic[site]
        site_pos = site_up - 1
        link_site_total_dic = {}

        E_value_list = []
        for i in range(site_up, site_down + 1):
            E_value = float(site_table[i].rstrip("\n").split(',')[8])
            E_value_list.append(E_value)
        Best_e_value = min(E_value_list)

        Pro1Name = get_linked_site_inform(site)[0]
        Pro2Name = get_linked_site_inform(site)[1]
        if Pro1Name == Pro2Name:
            Cross_type = "Intra"
        else:
            Cross_type = "inter"

        link_site_total_dic[site] = [
            site_correct(site), site_table[site_pos].strip("\n").split(',')[3],
            str(Best_e_value), site_table[site_up].rstrip("\n").split(',')[9],
            site_table[site_up].rstrip("\n").split(',')[5],
            site_table[site_up].rstrip("\n").split(',')[4], Cross_type
        ]

        raw_sub_spectra_dic = {}
        raw_sub_E_value_dic = {}
        raw_sub_peptide_dic = {}
        for raw in raw_name_list:
            raw_sub_spectra_dic[raw] = 0
            raw_sub_E_value_dic[raw] = []
            raw_sub_peptide_dic[raw] = []
            for j in range(site_up, site_down + 1):
                line_list = site_table[j].rstrip("\n").split(',')
                line_raw = line_list[2][:line_list[2].find(".")]
                if line_raw == raw:
                    raw_sub_spectra_dic[raw] += 1
                    raw_sub_E_value_dic[raw].append(float(line_list[8]))
                    raw_sub_peptide_dic[raw].append(line_list[5])
            if raw_sub_spectra_dic[raw] == 0:
                link_site_total_dic[site].append("")
                link_site_total_dic[site].append("")
                link_site_total_dic[site].append("")
            else:
                link_site_total_dic[site].append(str(raw_sub_spectra_dic[raw]))
                link_site_total_dic[site].append(
                    str(min(raw_sub_E_value_dic[raw])))
                link_site_total_dic[site].append(
                    str(len(list(set(raw_sub_peptide_dic[raw])))))
        if int(link_site_total_dic[site][1]) > spec_cutoff and float(
                link_site_total_dic[site][2]) < E_value_cutoff:
            b.write('\t'.join(link_site_total_dic[site]))
            b.write('\n')
    b.close()
    return


def statistic_report_file():
    c = open("report.txt", 'a')
    rep_table = open("report.txt").readlines()
    col_dic = {}
    total_spectra = 0
    total_colom = len(rep_table[0].strip().split("\t"))
    intra_num = 0
    for i in range(1, len(rep_table)):
        total_spectra += int(rep_table[i].strip("\n").split("\t")[1])
        if rep_table[i].strip("\n").split("\t")[5] == "intra":
            intra_num += 1

    col_dic[5] = float(intra_num) / (len(rep_table) - 1)
    col_dic[0] = len(rep_table) - 1
    col_dic[1] = total_spectra
    col_dic[2] = ""
    col_dic[3] = ""
    col_dic[4] = ""
    col_dic[6] = ""

    column_sub_dic = {}
    for k in [7, 9]:
        for j in range(k, total_colom, 3):
            column_sub_dic[j] = []
            for line in rep_table[1:]:
                line_list = line.rstrip("\n").split("\t")
                if line_list[j]:
                    column_sub_dic[j].append(int(line_list[j]))
                else:
                    continue
            col_dic[j] = sum(column_sub_dic[j])

    for j in range(8, total_colom, 3):
        column_sub_dic[j] = []
        for line in rep_table[1:]:
            line_list = line.rstrip("\n").split("\t")
            if line_list[j]:
                column_sub_dic[j].append(float(line_list[j]))
            else:
                continue
        col_dic[j] = str(min(column_sub_dic[j])) + 'To' + str(
            max(column_sub_dic[j]))
    last = []
    for k in range(total_colom):
        last.append(str(col_dic[k]))
    c.write("\t".join(last))
    print(last)
    c.close()
    return


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
        tab1 = open(link_site_file, 'r').readlines()
        contain_sameset = False
        for line in tab1:
            line_list_first_cell = line.rstrip("\n").split(",")[0]
            if line_list_first_cell == "SameSet":
                contain_sameset = True
                print("There are SameSets in this file")
                break
            else:
                continue
        if contain_sameset is False:
            site_pos_dic = get_site_lines(tab1)
            get_crosslink_site_info(tab1, site_pos_dic)
            statistic_report_file()
        else:
            remove_sameset(tab1)
            tmp_2 = open("tmpfile", 'r').readlines()
            site_pos_dic = get_site_lines(tmp_2)
            get_crosslink_site_info(tmp_2, site_pos_dic)
            statistic_report_file()
            os.remove("tmpfile")


if __name__ == "__main__":
    main()
