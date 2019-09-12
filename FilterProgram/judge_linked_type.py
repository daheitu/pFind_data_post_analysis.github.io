import os
import re
os.chdir(r"D:\E\Collabaration\SAGA complex\test")


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


f = open("list.txt", 'r').readlines()
b = open("list_type.txt", 'w')
for line in f:
    b.write("\t".join([line.strip(), judeg_linked_type(line.strip())]))
    b.write("\n")
    
b.close()
print("Done")