# coding = utf-8

import os, re
os.chdir(r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.05.08.15.12.32_AB_4D_trypsin\reports")
fltCerti = (3, 0.01)
trypAspn = r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.05.08.18.34.34_AB_RT_TRYP_ASPN\reports\ADLO_BSA_2019.05.08_Asp-N.Trypsin_KArGOv4.txt"
tryp = r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.05.08.17.05.35_AB_RT_trypsin\reports\ADLO_BSA_2019.05.08_Trypsin_KArGOv4.txt"
#trypAspn = r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.05.08.22.39.05_AB_M20_TRY_ASP\reports\ADLO_BSA_2019.05.08_Asp-N.Trypsin_KArGOv4.txt"
#tryp = r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.05.08.20.38.03\reports\ADLO_BSA_2019.05.08_Trypsin_KArGOv4.txt"
#trypAspn = r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.05.08.16.09.51_AB_4D_trypsin_ASPN\reports\ADLO_BSA_2019.05.08_Asp-N.Trypsin_KArGOv4.txt"
#tryp = r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.05.08.15.12.32_AB_4D_trypsin\reports\ADLO_BSA_2019.05.08_Trypsin_KArGOv4.txt"

#f = open("ADLO_BSA_2019.05.08_Trypsin_KArGOv4.txt", 'r').readlines()


def judgeProtType(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]

    if protein1 != protein2:
        return "Inter"
    else:
        return "Intra"


def mergeCond(filePath):
    f = open(filePath, 'r').readlines()
    condInfoDic = {}
    columnList = f[0].strip().split("\t")
    for k in range(6, len(columnList),3):
        print(k)
        rawName = "_".join(columnList[k].split("_")[:-1])
        condTion = "_".join(columnList[k].split("_")[:4])
        print(condTion)
        if condTion not in condInfoDic:
            condInfoDic[condTion] = {}

            for line in f[1:]:
                lineList = line.rstrip("\n").split("\t")
                linkPair = lineList[0]
                protType = judgeProtType(linkPair)
                if protType != "Intra":
                    continue
                else:
                    if lineList[k] != "":
                        spec = int(lineList[k])
                        evalue = float(lineList[k+1])
                        condInfoDic[condTion][linkPair] = [spec, evalue]
        else:
            for line in f[1:]:
                lineList = line.rstrip("\n").split("\t")
                linkPair = lineList[0]
                if lineList[k] != "":
                    spec = int(lineList[k])
                    evalue = float(lineList[k+1])
                    if linkPair not in condInfoDic[condTion]:
                        condInfoDic[condTion][linkPair] = [spec, evalue]
                    else:
                        condInfoDic[condTion][linkPair][0] += spec
                        condInfoDic[condTion][linkPair][1] = min(condInfoDic[condTion][linkPair][1], evalue)
    return condInfoDic


def mergeTryAsp(tryDic, aspDic):
    mergeDic = {}
    for key in tryDic:
        mergeDic[key] = {}
        allLinkPair = list(set(tryDic[key].keys()) | set(aspDic[key].keys()))
        allLinkPair.sort()
        print(allLinkPair)
        for pair in allLinkPair:
            if pair in tryDic[key] and pair in aspDic[key]:
                spec =  tryDic[key][pair][0] + aspDic[key][pair][0]
                evalue = min(tryDic[key][pair][1], aspDic[key][pair][1])
                mergeDic[key][pair] = [spec, evalue]
            elif pair in tryDic[key]:
                mergeDic[key][pair] = tryDic[key][pair]
            else:
                mergeDic[key][pair] = aspDic[key][pair]
    return mergeDic


def filterDic(tryDic, certificate):
    specCutOFF, evalueCutOFF = certificate
    fltDic = {}
    for cond in tryDic:
        fltDic[cond] = {}
        for pair in tryDic[cond]:
            spec, evalue = tryDic[cond][pair]
            if spec >= specCutOFF and evalue <= evalueCutOFF:
                fltDic[cond][pair] = tryDic[cond][pair]
    
    return fltDic
    

def calNumOfLinkSite(repDic):
    for cond in repDic:
        numOfXS = len(repDic[cond])
        numOfSpec = 0
        for pair in repDic[cond]:
            numOfSpec += repDic[cond][pair][0]
        print(cond, numOfXS, numOfSpec)


def main():
    trypCondInfoDic = mergeCond(tryp)
    trypAspnCondInfoDic = mergeCond(trypAspn)
    trypFltDic = filterDic(trypCondInfoDic, fltCerti)
    trypAspnFltDic = filterDic(trypAspnCondInfoDic, fltCerti)
    mergeDic = mergeTryAsp(trypFltDic, trypAspnFltDic)
    
    calNumOfLinkSite(mergeDic)
    #mergeDic = mergeTryAsp(trypCondInfoDic, trypAspnCondInfoDic)
    #fltDic = filterDic(mergeDic, (3, 0.01))
    #calNumOfLinkSite(fltDic)


if __name__ == "__main__":
    main()
