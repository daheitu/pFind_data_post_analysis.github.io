import os
import re
import csv


os.chdir(r"C:\Users\Yong\Desktop\xiNET")

delta_dic = {"Utp8 ": 300, "Utp4": 0, "Utp15 ": 380, "Utp9 ": 360, "Utp5 " : 432}
# dict() csv.reader()


def site_correct(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]
    return protein1, protein2, position1, position2


rep = open("kargo_utpa_correct.txt", 'w')
with open("kargo_utpa.txt", 'r') as cscfile:
    f = csv.reader(cscfile, delimiter='\t', quotechar='\"')
    for row in f:        
        print(row[0])
        site = row[0]
        prot1, prot2, pos1, pos2 = site_correct(site)
        pos1 = int(pos1) + delta_dic[prot1]
        pos2 = int(pos2) + delta_dic[prot2]

        new_site = prot1+"("+ str(pos1) + ")-" + prot2 + "(" + str(pos2) + ")"
        rep.write("\t".join([row[0], new_site]))
        rep.write("\n")

rep.close()
