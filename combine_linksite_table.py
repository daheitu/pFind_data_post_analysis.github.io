import os
import re
os.chdir(r"D:\E\Collabaration\SheYang\results")
file_list = os.listdir(os.getcwd())


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
        if abs(int(position1)-int(position2)) < 5:
            link_type = "Inter"
        else:
            link_type = "Intra" 
    return protein1, protein2, position1, position2, link_type


def judeg_linked_type(linked_site):
    if "/" in linked_site:
        site_list = linked_site.split("/")
        for site in site_list:
            if "REVERSE" in site:
                site_list.remove(site)
            else:
                continue
        
        type_list = []
        for site in site_list:
            type_list.append(get_linked_site_inform(site)[-1])
        
        if "Intra" in type_list:
            return "Intra"
        else:
            return "Inter"
    else:
        return get_linked_site_inform(linked_site)[-1]


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


def site_list_correction(linked_site):
    if "/" in linked_site:
        site_list = linked_site.split("/")
        for site in site_list:
            if "REVERSE" in site:
                site_list.remove(site)
            else:
                continue
        for i in range(len(site_list)):
            site_list[i] = site_correct(site_list[i])
        site_list.sort()
        return "/".join(site_list)
    else:
        return site_correct(linked_site)


file_site_dic = {}

file_list = os.listdir(os.getcwd())

for fl in file_list:
    site_dic = {}
    f = open(fl, 'r').readlines()
    for line in f[1:-1]:
        if line == "":
            continue
        else:
            line_list = line.rstrip("\n").split("\t")
            site = site_list_correction(line_list[0])

            if site not in site_dic:
                site_dic[site] = line_list[1:3]
            else:
                print(line)
                print("wrong")
    file_site_dic[fl] = site_dic

site_list = []
for fl in file_site_dic:
    for site in file_site_dic[fl]:
        if site not in site_list:
            site_list.append(site)


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
