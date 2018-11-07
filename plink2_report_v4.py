# coding: utf-8
"""
Created on Mon Aug 15 10:38:11 2017
@author: CaoYong
This script can help you to summary the plink2 report file
"""

import os
import re

os.chdir(
    r"C:\Users\Yong\Documents\pLink\pLink_task_2018.11.07.10.52.00\reports"
)
spec_cutoff = 0  # spectra number cut-off
Best_evalue_cutoff = 2
E_value_cutoff_SpecLvl = 1


def site_correct(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]

    if protein1 != protein2:
        link_type = "Inter"
    else:
        if abs(int(position1) - int(position2)) < 5:
            link_type = "Inter"
        else:
            link_type = "Intra"

    if int(position1) <= int(position2):
        correc_site = linked_site
    else:
        a = linked_site.split(")-")[0] + ")"
        b = linked_site.split(")-")[1]
        correc_site = b + "-" + a

    return correc_site, link_type


def site_list_process(site_list):
    i = 0
    while i < len(site_list):
        if "REVERSE" in site_list[i] or "gi|CON" in site_list[i]:
            site_list.remove(site_list[i])
        else:
            i += 1

    link_type_list = []
    for i in range(len(site_list)):
        site_list[i] = site_correct(site_list[i])[0]
        link_type = site_correct(site_list[i])[1]
        if link_type not in link_type_list:
            link_type_list.append(link_type)
        else:
            pass
    site = "/".join(site_list)
    if link_type_list == ["Inter"]:
        return site, "Inter"
    else:
        return site, "Intra"


def get_report_file_name():
    path = os.getcwd()
    path_d = os.path.dirname(path)
    os.chdir(path_d)
    print(path_d)
    file_list = os.listdir(path_d)
    for fl in file_list:
        if fl[-5:] in ["pfind", "plink"]:
            para = open(fl).readlines()
            for line in para:
                if line[:10] == "spec_title":
                    print(line)
                    spec_title = line.rstrip("\n").split("=")[1].strip()
                elif line[:11] == "enzyme_name":
                    enzyme = line.rstrip("\n").split("=")[1].strip()
                elif line[:7] == "linker1":
                    linker = line.rstrip("\n").split("=")[1].strip()
                else:
                    continue
        else:
            continue

    report_file_name = spec_title + "_" + enzyme + "_" + linker + "v4.txt"
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
        "Peptide", "Prote type"
    ]
    for name in raw_name_list:
        col.append(name + "_SpecNum")
        col.append(name + "_E-value")
        col.append(name + "_UniquePepNum")

    b.write('\t'.join(col))
    b.write('\n')
    return raw_name_list


def statistic_report_file():
    report_file_name = get_report_file_name()
    c = open(report_file_name, 'a')
    rep_table = open(report_file_name, 'r').readlines()
    title_list = rep_table[0].split("\t")
    col_dic = {}
    total_spectra = 0
    total_colom = len(rep_table[0].strip().split("\t"))
    intra_num = 0
    for i in range(1, len(rep_table)):
        total_spectra += int(rep_table[i].strip("\n").split("\t")[1])
        if rep_table[i].strip("\n").split("\t")[5] == "Intra":
            intra_num += 1
    print(intra_num)
    col_dic[5] = float(intra_num) / (float(len(rep_table)) - 1.0)
    col_dic[0] = len(rep_table) - 1
    col_dic[1] = total_spectra
    col_dic[2] = ""
    col_dic[3] = ""
    col_dic[4] = ""

    column_sub_dic = {}
    for k in [6, 8]:
        for j in range(k, total_colom, 3):
            col_dic[j] = 0
            for line in rep_table[1:]:
                line_list = line.rstrip("\n").split("\t")
                if line_list[j]:
                    col_dic[j] += int(line_list[j])
                else:
                    continue

    for j in range(7, total_colom, 3):
        line_1_list = rep_table[1].rstrip("\n").split("\t")
        column_sub_dic[j] = []
        min_evl = float(line_1_list[j])
        max_evl = float(line_1_list[j])
        for line in rep_table[2:]:
            line_list = line.rstrip("\n").split("\t")
            evalue = line_list[j]
            if evalue:
                evalue = float(line_list[j])
                if evalue > max_evl:
                    max_evl = evalue
                elif evalue < min_evl:
                    min_evl = evalue
                else:
                    pass
            else:
                continue
        col_dic[j] = str(min_evl) + ' To ' + str(max_evl)
    last = []
    for k in range(total_colom):
        last.append(str(col_dic[k]))
    c.write("\t".join(last))
    c.write("\n")
    print(last)
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
        if f[p].rstrip("\n").split(",")[0] in ["SameSet", "SubSet"]:
            while p < len(f) and f[p].rstrip("\n").split(",")[0] in [
                    "SameSet", "SubSet"
            ]:
                if f[p].rstrip("\n").split(",")[0] == "SameSet":
                    site_list.append(f[p].rstrip("\n").split(",")[1])
                else:
                    pass
                p += 1
        else:
            pass
        site = site_list_process(site_list)[0]
        if site == "":
            while p < len(f) and f[p].rstrip("\n").split(",")[0] == "":
                p += 1
        else:
            link_type = site_list_process(site_list)[1]
            best_svm_score = f[p].rstrip("\n").split(",")[9]
            pep_std = f[p].rstrip("\n").split(",")[5]
            while p < len(f) and f[p].rstrip("\n").split(",")[0] == "":
                line_list = f[p].rstrip("\n").split(",")
                evalue = line_list[8]
                if float(evalue) < E_value_cutoff_SpecLvl:
                    raw_name = line_list[2][:line_list[2].find(".")]
                    pep = line_list[5]
                    if raw_name not in spec_dic:
                        spec_dic[raw_name] = 1
                    else:
                        spec_dic[raw_name] += 1

                    if raw_name not in pep_dic:
                        pep_dic[raw_name] = [pep]
                    else:
                        if pep not in pep_dic[raw_name]:
                            pep_dic[raw_name].append(pep)
                        else:
                            pass

                    if raw_name not in evalue_dic:
                        evalue_dic[raw_name] = float(evalue)
                    else:
                        if float(evalue) < evalue_dic[raw_name]:
                            evalue_dic[raw_name] = float(evalue)
                        else:
                            pass
                else:
                    pass
                p += 1

            total_spec = 0
            min_evalue = 1
            for key in evalue_dic:
                total_spec += spec_dic[key]
                if evalue_dic[key] < min_evalue:
                    min_evalue = evalue_dic[key]
                else:
                    continue

            rep_dic = {}
            for raw_name in raw_name_list:
                if raw_name not in spec_dic:
                    rep_dic[raw_name] = ["", "", ""]
                else:
                    rep_dic[raw_name] = [
                        str(spec_dic[raw_name]),
                        str(evalue_dic[raw_name]),
                        str(len(pep_dic[raw_name]))
                    ]
            rep_list = [
                site, str(total_spec),
                str(min_evalue), best_svm_score, pep_std, link_type
            ]
            for raw_name in raw_name_list:
                rep_list.extend(rep_dic[raw_name])

            if float(rep_list[1]) > spec_cutoff and float(
                    rep_list[2]) < Best_evalue_cutoff:
                b.write("\t".join(rep_list))
                b.write("\n")
            else:
                pass

        n = p
    b.close()
    # statistic_report_file()


if __name__ == "__main__":
    main()
