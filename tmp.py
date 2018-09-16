import os
import re

path = r"C:\Users\Yong\Desktop\DISTANCE\Cal_Pdb_Distance\CNGP"
XL_sites_list = ["R", "K"]
os.chdir(path)


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
        if abs(int(position1) - int(position2)) < 2:
            link_type = "Inter"
        else:
            link_type = "Intra"

    if int(position1) <= int(position2):
        correc_site = linked_site
    else:
        a = linked_site.split(")-")[0] + ")"
        b = linked_site.split(")-")[1]
        correc_site = b + "-" + a

    return link_type

c = open("judege_type.txt", 'w')
f = open("list.txt").readlines()
pair_list = []
for line in f[1:]:
    line_list = line.strip().split("\t")
    if line_list[0].isdigit():
        break
    else:
        pair_list.append(line_list[0])
for pairs in pair_list:
    c.write("\t".join([pairs, site_correct(pairs)]))
    c.write("\n")

c.close()
