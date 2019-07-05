# coding = utf8
'''
python 3.6, pLink2.3.5 or later version
if you have any question, please contact caoyong@nibs.ac.cn
'''
import os
import re

scoreType = 0 # important, 0 for spectra num, 1 for Best E-value, 2 for Best_svm
path = r'D:\softwareData\plink2\pLink_task_2018.03.05.07.48.38_KARGO_T\reports' # the path of reports file

###############Don't change the following lines#####################################

os.chdir(path)


def splitLinksite(linked_site):
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
        if abs(int(position1)-int(position2)) < 3:
            link_type = "Inter"
        else:
            link_type = "Intra" 
    return protein1, protein2, position1, position2


def getCrossLinkSiteInfo(openfl):
    siteInfoDic = {}
    i = 2
    while i < len(openfl):
        lineList = openfl[i].rstrip('\n').split(",")
        if lineList[0].isdigit():
            sitePair = lineList[1]
            totalSpecNum = lineList[-1]
        else:
            print('wrong')
        
        svmScoreList = []
        eValueList = []
        p = i + 1
        while p < len(openfl):
            lineList = openfl[p].rstrip('\n').split(",")
            if lineList[0] in ["SameSet", "SubSet"]:
                p += 1
            elif lineList[0].isdigit():
                break
            else:
                eValueList.append(float(lineList[8]))
                svmScoreList.append(float(lineList[9]))
                p += 1
        bestSvmScore = min(svmScoreList)
        bestEvalue = min(eValueList)
        prot1, prot2, pos1, pos2 = splitLinksite(sitePair)
        if scoreType == 0:
            siteInfoDic[sitePair] = [totalSpecNum, prot1, prot2, pos1, pos2]
        elif scoreType == 1:
            siteInfoDic[sitePair] = [bestEvalue, prot1, prot2, pos1, pos2]
        elif scoreType == 2:
            siteInfoDic[sitePair] = [bestSvmScore, prot1, prot2, pos1, pos2]
        else:
            print("scoreType is wrong")
        
        i = p
    return siteInfoDic
        

def main():
    for fl in os.listdir(os.getcwd()):
        if fl[-22:] == "cross-linked_sites.csv":
            f = open(fl).readlines()
            rep = open(fl[:-22] + "_xiNET.csv", 'w')
            rep.write(",".join(["Score", "Protein1", "Protein2",\
                    "LinkPos1", "LinkPos2"])+'\n')
            siteInfoDic = getCrossLinkSiteInfo(f)
            for site in siteInfoDic:
                rep.write(",".join(siteInfoDic[site])+'\n')
            rep.close()
        else:
            continue


if __name__ == "__main__":
    main()