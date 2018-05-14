import os
import re
os.chdir(r"C:\Users\Yong\Desktop\results\combine")


def gene_blank_list(list_length):
    blank_list = []
    for i in range(list_length):
        blank_list.append("")
    return blank_list


def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find("-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:n].strip()
    if protein1 != protein2:
        link_type = "Inter"
    else:
        if abs(int(position1) - int(position2)) < 5:
            link_type = "Inter"
        else:
            link_type = "Intra"
    return protein1, protein2, position1, position2, link_type


def judeg_linked_type(linked_site):
    if "/" in linked_site:
        site_list = linked_site.split("/")
        type_list = []
        for site in site_list:
            type_list.append(get_linked_site_inform(site)[-1])

        if "Intra" in type_list:
            return "Intra"
        else:
            return "Inter"
    else:
        return get_linked_site_inform(linked_site)[-1]


def combine_table(path):
    file_list = os.listdir(path)
    file_site_dic = {}
    file_line_length_dic = {}

    for fl in file_list:
        site_dic = {}
        f = open(fl, 'r').readlines()
        line_length = len(f[1].rstrip("\n").split("\t")) - 1
        for line in f:
            if line == "":
                continue
            else:
                line_list = line.rstrip("\n").split("\t")
                site = line_list[0]
                if site not in site_dic:
                    site_dic[site] = line_list[1:]
                else:
                    print(line)
                    print("wrong")
        file_site_dic[fl] = site_dic
        file_line_length_dic[fl] = line_length

    full_site_list = []
    for fl in file_site_dic:
        for site in file_site_dic[fl]:
            if site not in full_site_list:
                full_site_list.append(site)
            else:
                continue

    table_dic = {}
    for site in full_site_list:
        table_dic[site] = [site]
        for fl in file_site_dic:
            if site in file_site_dic[fl]:
                table_dic[site].extend(file_site_dic[fl][site])
            else:
                table_dic[site].extend(
                    gene_blank_list(file_line_length_dic[fl]))
    return table_dic


def main():
    table_dic = combine_table(os.getcwd())
    rep = open("report.txt", "w")
    for site in table_dic:
        table_dic[site].append(judeg_linked_type(site))
        rep.write("\t".join(table_dic[site]))
        rep.write("\n")

    rep.close()
    print("done")


if __name__ == "__main__":
    main()

"""
rep = open("report.txt", "w")
title_list = ["Site Paire"]
table_dic = {}
for fl in file_site_dic:
        title_list.extend([fl+"_SpecNum", fl + "_E-value"])
title_list.append("Linked_type")
rep.write("\t".join(title_list))
rep.write("\n")
for site in site_list:
    table_dic[site] = [site]
    for fl in file_site_dic:
        if site in file_site_dic[fl]:
            table_dic[site].extend(file_site_dic[fl][site])
        else:
            table_dic[site].extend(["", ""])
    table_dic[site].append(judeg_linked_type(site))
    
    rep.write("\t".join(table_dic[site]))
    rep.write("\n")
rep.close()
print("done")
"""