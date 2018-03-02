import os
import re

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


os.chdir(r"D:\E\Collabaration\SAGA complex\test")
f = open("report.txt", 'r').readlines()
b = open("SAGA xiNET.csv", "w")
b.write(",".join(["Score", "Protein1", "Protein2", "LinkPos1", "LinkPos2"]))
b.write("\n")


for line in f[1:]:
    line_list = line.rstrip("\n").split("\t")
    site = line_list[0]
    spectra = line_list[1]
    write_list = [spectra]
    if "/" in site:
        continue
    else:
        if "Molecular" in site:
            continue
        else:
            write_list.extend(get_linked_site_inform(site)[:-1])
            b.write(",".join(write_list))
            b.write("\n")
b.close()
print("Done")
