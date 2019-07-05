# coding = utf-8

import os, re
#import plotly.plotly as py
#import plotly.graph_objs as go
 

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
#c = open("ALPK1_scatterbuble_N.txt", 'w')
b.write("\t".join(["label","x_aisx", "y_asix", "spec", "Link_type"])+ "\n")
#c.write("\t".join(["x_aisx", "y_asix", "spec"])+ "\n")
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
        if int(pos1) < int(pos2):
            wList2 = [str(f.index(line)), pos1, pos2, spec_N, "N"]  
            wList1 = [str(f.index(line))+ "'", pos2, pos1, spec_NC, "NC"]
        else:
            wList2 = [str(f.index(line)), pos2, pos1, spec_N, "N"]  
            wList1 = [str(f.index(line))+ "'", pos1, pos2, spec_NC, "NC"]
        if "0.01" not in wList1:
            b.write("\t".join(wList1) + "\n")
        if "0.01" not in wList2:
            b.write("\t".join(wList2) + "\n")
        
b.close()

def import_data():
    x_axis = []
    y_axis = []
    spec_list = []
    label_list = []
    f = open("ALPK1_scatterbuble.txt", 'r').readlines()
    
    for line in f[1:]:
        lineList = line.strip().split("\t")
        x_axis.append(int(lineList[0]))
        y_axis.append(int(lineList[1]))
        spec_list.append(int(lineList[2]))
        label_list.append(f.index(line))
    return x_axis, y_axis, spec_list, label_list

