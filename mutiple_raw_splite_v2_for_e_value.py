"""
Created on Mon Aug 15 10:38:11 2016
@author: CaoYong
"""


import os
os.chdir(
    r"D:\softwareData\plink2\20160703\search_task_2016.07.05.15.18.33_Aldo_3\reports"
)


# 对list去重复
def reduce_redundancy_list(list):
    new_site_list = []
    for word in list:
        if word not in new_site_list:
            new_site_list.append(word)
    return new_site_list



# 找到包含site信息的行
def get_site_lines(filename):
    site_table = open(filename).readlines()
    ordr_pos = []
    for i in range(len(site_table)):
        if len([
                x for x in site_table[i].strip().split(',') if x != ''
        ]) == 4:
            ordr_pos.append(i)
    ordr_pos.append(len(site_table))
    return ordr_pos, site_table


# 去掉sameset
def remove_sameset():
    filename = os.listdir(os.getcwd())
    for fl in filename:
        if fl[-22:-4] == "cross-linked_sites":
            ordr_pos, site_table = get_site_lines(fl)
    tmp = open("tmpfile", "w")
    for i in range(1,len(ordr_pos) - 1):
        Con_sameset = False
        list_sameset = []
        list_sameset.append(site_table[ordr_pos[i]].rstrip("\n").split(",")[1])
        for j in range(ordr_pos[i], ordr_pos[i + 1]):
            first_cell = site_table[j].rstrip("\n").split(",")[0]
            if first_cell == "Sameset":
                Con_sameset = True
        if Con_sameset == False:
            for j in range(ordr_pos[i], ordr_pos[i + 1]):
                tmp.write(site_table[j])
        else:
            for j in range(ordr_pos[i], ordr_pos[i + 1]):
                first_cell = site_table[j].rstrip("\n").split(",")[0]
                if first_cell == "Sameset":
                    list_sameset.append(
                        site_table[j].rstrip("\n").split(",")[1])
                elif first_cell == "":
                    break
            link_site = "/".join(list_sameset)
            site_info_line_list = site_table[ordr_pos[i]].rstrip("\n").split(
                ",")
            site_info_line_list[1] = link_site
            tmp.write(",".join(site_info_line_list))
            tmp.write("\n")
            for j in range(ordr_pos[i] + 1, ordr_pos[i + 1]):
                first_cell = site_table[j].rstrip("\n").split(",")[0]
                if first_cell != "Sameset":
                    tmp.write(site_table[j])
    tmp.close()
    return


def get_crosslink_site_info():
    ordr_pos, site_table_new = get_site_lines("tmpfile")
    b = open("report.txt", 'w')
    sit_list = []
    raw_name_list = []
    for m in range(len(ordr_pos) - 1):
        for k in range(ordr_pos[m] + 1, ordr_pos[m + 1]):
            a = site_table_new[k].strip().split(',')[2].find(".")
            raw_name_list.append(site_table_new[k].strip().split(',')[2][:a])
            raw_name_list = reduce_redundancy_list(raw_name_list)
    raw_name_list.sort()
    print(raw_name_list)
    col = [
        "Linked_site", "Total_spec", "Best_e_value", "Best_svm", "Peptide",
        "Peptide_mass", "Prote_type"
    ]
    for name in raw_name_list:
        col.append(name + "_spec")
        col.append(name + "_E_value")
        col.append(name + "unique_pep_num")

    b.write(','.join(col))
    b.write('\n')
    m = 0
    link_site_total_dic = {}
    for l in range(len(ordr_pos) - 1):
        linked_site = site_table_new[ordr_pos[l]].strip().split(',')[1]
        E_value_list = []
        for i in range(ordr_pos[l] + 1, ordr_pos[l + 1]):
            E_value_list.append(float(site_table_new[i].strip().split(',')[8]))
        E_value_list = reduce_redundancy_list(E_value_list)
        if "0" in E_value_list:
            E_value_list.remove(0)
        Best_e_value = min(E_value_list)
        link_site_total_dic[linked_site] = [
            site_table_new[ordr_pos[l]].strip().split(',')[1],
            site_table_new[ordr_pos[l]].strip().split(',')[3],
            str(Best_e_value),
            site_table_new[ordr_pos[l] + 1].strip().split(',')[9],
            site_table_new[ordr_pos[l] + 1].strip().split(',')[5],
            site_table_new[ordr_pos[l] + 1].strip().split(',')[4]
        ]
        m = linked_site.find("(")
        n = linked_site.find(")")
        p = linked_site.find("-")
        x = linked_site.find("(", p)
        y = linked_site.find(")", p)
        if linked_site[:m] == linked_site[p + 1:x]:
            link_site_total_dic[linked_site].append("intra")
        else:
            link_site_total_dic[linked_site].append("inter")
        raw_sub_spectra_dic = {}
        raw_sub_E_value_dic = {}
        raw_sub_peptide_dic = {}
        for raw in raw_name_list:
            raw_sub_spectra_dic[raw] = 0
            raw_sub_E_value_dic[raw] = []
            raw_sub_peptide_dic[raw] = []
            for j in range(ordr_pos[l] + 1, ordr_pos[l + 1]):
                if raw == site_table_new[j].strip().split(',')[
                        2][:site_table_new[j].strip().split(',')[2].find(".")]:
                    raw_sub_spectra_dic[raw] += 1
                    raw_sub_E_value_dic[raw].append(
                        float(site_table_new[j].strip().split(',')[8]))
                    raw_sub_peptide_dic[raw].append(
                        site_table_new[j].strip().split(',')[5])
            if raw_sub_spectra_dic[raw] == 0:
                link_site_total_dic[linked_site].append("")
                link_site_total_dic[linked_site].append("")
                link_site_total_dic[linked_site].append("")
            else:
                link_site_total_dic[linked_site].append(
                    str(raw_sub_spectra_dic[raw]))
                link_site_total_dic[linked_site].append(
                    str(min(raw_sub_E_value_dic[raw])))
                link_site_total_dic[linked_site].append(
                    str(len(list(set(raw_sub_peptide_dic[raw])))))
        b.write(','.join(link_site_total_dic[linked_site]))
        b.write('\n')
    b.close()
    return


def statistic_report_file():
    c = open("report.txt", 'a')
    rep_table = open("report.txt").readlines()
    col_dic = {}
    total_spectra = 0
    total_colom = len(rep_table[0].strip().split(","))
    intra_num = 0
    for i in range(1, len(rep_table)):
        total_spectra += int(rep_table[i].strip("\n").split(",")[1])
        if rep_table[i].strip("\n").split(",")[5] == "intra":
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
            for i in range(1, len(rep_table)):
                if rep_table[i].strip("\n").split(",")[j]:
                    column_sub_dic[j].append(
                        int(rep_table[i].strip("\n").split(",")[j]))
            col_dic[j] = sum(column_sub_dic[j])

    for j in range(8, total_colom, 3):
        column_sub_dic[j] = []
        for i in range(1, len(rep_table)):
            if rep_table[i].strip("\n").split(",")[j]:
                column_sub_dic[j].append(
                    round(float(rep_table[i].strip("\n").split(",")[j]), 2))
        col_dic[j] = str(
            str(min(column_sub_dic[j])) + ',' + str(max(column_sub_dic[j])))
    last = []
    for k in range(total_colom):
        last.append(str(col_dic[k]))
    c.write(",".join(last))
    print(last)
    c.close()
    return


def main():
    remove_sameset()
    get_crosslink_site_info()
    statistic_report_file()
    os.remove("tmpfile")


if __name__ == "__main__":
    main()
