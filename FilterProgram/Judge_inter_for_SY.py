import os
import re
os.chdir(r"C:\Users\Yong\Desktop\sgc")


def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    # p = linked_site.find("-")
    # m = linked_site.find("(" + position1 + ")-")
    # n = linked_site.find("(" + position2 + ")", p)
    # protein1 = linked_site[:m].strip()
    # protein2 = linked_site[p + 1:n].strip()
    if int(position1) < 500 and int(position2) > 900:
        return "Inter"
    else:
        return "Intra"


f = open("Sheyang list.txt", 'r').readlines()
b = open("report.txt", 'w')
for line in f:
    site = line.strip()
    link_type = get_linked_site_inform(site)
    b.write("\t".join([site, link_type]))
    b.write("\n")

b.close()
