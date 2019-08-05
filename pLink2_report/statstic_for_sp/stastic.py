# coding = utf8

import os, re
 
os.chdir(r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.07.03.15.20.56\reports")

def getSiteInfo(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]
    return protein1, protein2, position1, position2


def statsticColumn(k, rawName, openedFL, flToWrite):
    f = openedFL
    b = flToWrite
    numOfPepOfNTxl = 0
    numOfSpecOfNTxl = 0
    numOfPepOfMMxl = 0
    numOfSpecOfMMxl = 0
    numOfPepfreeOfMMxl = 0
    numOfSpecfreeOfMMxl = 0
    numOfPepAcIR7Oxl = 0
    numOfSpecAcIR7Oxl = 0
    numOfPepEL10Oxl = 0
    numOfSpecEL10xl = 0
    numOfPepOthOxl = 0
    numOfSpecOthxl = 0
    for line in f[1:]:
        #print(line.strip().split("\t"))
        lineList = line.rstrip("\n").split("\t")
        linkSite = lineList[0]
        pro1, pro2, pos1, pos2 = getSiteInfo(linkSite)
        if lineList[k] not in ["", "1"]:#]:
            totalSpec = int(lineList[k])
            if pro1[-4:] == "PR-4" and pro2[-4:] == "PR-4":
                if pos1 == "1" or pos2 == "1":
                    numOfPepOfNTxl += 1
                    numOfSpecOfNTxl += totalSpec
                else:
                    numOfPepOfMMxl += 1
                    numOfSpecOfMMxl += totalSpec
            elif pro1 == pro2:
                numOfPepfreeOfMMxl += 1
                numOfSpecfreeOfMMxl += totalSpec
            else:
                if "AcIR7(5)" in linkSite and "(4)" in linkSite:
                    numOfPepAcIR7Oxl += 1
                    numOfSpecAcIR7Oxl += totalSpec
                elif "EL10" in linkSite and "(4)" in linkSite:
                    numOfPepEL10Oxl += 1
                    numOfSpecEL10xl += totalSpec
                else:
                    numOfPepOthOxl += 1
                    numOfSpecOthxl += totalSpec
    wlist = [rawName, numOfPepOfNTxl, numOfSpecOfNTxl, numOfPepOfMMxl, \
            numOfSpecOfMMxl, numOfPepfreeOfMMxl, numOfSpecfreeOfMMxl, \
            numOfPepAcIR7Oxl, numOfSpecAcIR7Oxl, numOfPepEL10Oxl,\
            numOfSpecEL10xl, numOfPepOthOxl, numOfSpecOthxl]    
    b.write("\t".join([str(ele) for ele in wlist]) + "\n")


def main():
    f = open("AR7_con_2019.07.03_Trypsin_DSSv4.txt", 'r').readlines()
    b = open("report.txt", 'w')
    b.write("\t".join(["rawName", "numOfPepOfNTxl", "numOfSpecOfNTxl",\
        "numOfPepOfMMxl", "numOfSpecOfMMxl", "numOfPepfreeOfMMxl", \
        "numOfSpecfreeOfMMxl", "AcIR7O_Inter_sites", "AcIR7O_Inter_Spec",\
        "EL10_Inter_sites","EL10_Inter_Spec", "other_Inter_sites",\
        "other_Inter_Spec"])+"\n")
    titleList = f[0].strip().split("\t")
    for k in range(6, len(titleList), 3):
        rawName = titleList[k]
        statsticColumn(k, rawName, f, b)
    b.close()

if __name__ == "__main__":
    main()