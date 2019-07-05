# coding = utf-8

import os, re
 

os.chdir(r"E:\Script\scaterpiePlot\preTreatment")

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

f = open("scatterPie.csv", 'r').readlines()
b = open("ALPK1_scatterbuble.txt", 'w')
b.write("\t".join(["x_aisx", "y_asix", "spec"])+ "\n")
for line in f[1:]:
    lineList = line.rstrip("\n").split(",")
    print(line.rstrip("\n").split(","))
    linkPair = lineList[0]
    if linkPair == "":
        continue
    else:
        pos1, pos2 = site_correct(linkPair)[-2:]
        spec_NC = lineList[1]
        spec_N = lineList[3]
        wList1 = [pos1, pos2, spec_NC]
        b.write("\t".join(wList1) + "\n")
        wList2 = [pos2, pos1, spec_N]
        b.write("\t".join(wList2) + "\n")
        
b.close()