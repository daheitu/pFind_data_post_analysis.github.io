import os

os.chdir(r"G:\6peptide_PhGO_HF_20190223")

f1 = open("DK10_C_highBpep_list.txt", 'r').readlines()
f2 = open("DK10_PHGO_highBpep_list.txt", 'r').readlines()



def generate_mass_range(num, delta_ppm):
    delta = num * delta_ppm / 1000000
    return num - delta, num + delta


def getFileInfo(f):
    tableDic = {}
    for line in f[1:]:
        lineList = line.strip().split("\t")
        mz = float(lineList[0])
        tableDic[mz] = lineList
    
    return tableDic


def findOneMZ1(mz, ms1MzList,  tol_ppm):
    k = 0
    [low_ms, up_ms] = generate_mass_range(mz, tol_ppm)
    while k < len(ms1MzList):
        if float(ms1MzList[k]) > up_ms:
            return False, -1
        elif float(ms1MzList[k]) >= low_ms:
            return True, ms1MzList[k]
        else:
            k += 1
    if k == len(ms1MzList):
        return False, -1


def combineTwoTable(info1Dic, info2Dic):
    finalDic = {}
    for mz in info1Dic:
        if mz:
            info1BlankList = [' '] * len(info1Dic[mz])
            break
    for mz in info2Dic:
        if mz:
            info2BlankList = [' '] * len(info2Dic[mz])
            break
    print(info1BlankList, info2BlankList)
    mz1List = sorted(list(info1Dic.keys()))
    mz2List = sorted(list(info2Dic.keys()))
    #print(mz1List)

    for mz in mz1List:
        findBool, matchMZ = findOneMZ1(mz, mz2List, 10)
        if findBool:
            areaRatio = float(info2Dic[matchMZ][-1])/float(info1Dic[mz][-1])
            maxInts = max(float(info1Dic[mz][-3]), float(info2Dic[matchMZ][-3]))
            finalDic[mz] = info1Dic[mz] + info2Dic[matchMZ] + [round(areaRatio, 3)] + [maxInts]
            #mz2List.remove(matchMZ)
        else:
            finalDic[mz] = info1Dic[mz]+ info2BlankList
    
    for mz in mz2List:
        findBool, matchMZ = findOneMZ1(mz, mz1List, 10)
        if findBool:
            pass
        else:
            deltaArea = info2Dic[mz][-1]
            finalDic[mz] = info1BlankList + info2Dic[mz] + [deltaArea] + [info2Dic[mz][-3]]
    
    return finalDic


def writeDicToFile(wdic, flToW):
    mzList = sorted(wdic.keys())
    for key in mzList:
        wList = [str(ele) for ele in wdic[key]]
        flToW.write(",".join(wList) + "\n")


def main():
    f1dic = getFileInfo(f1)
    f2dic = getFileInfo(f2)
    #print(f1dic)
    combDic = combineTwoTable(f1dic, f2dic)
    #print(combDic)
    b = open("DK10_report_combine.csv", 'w')
    writeDicToFile(combDic, b)
    b.close()


if __name__ == "__main__":
    main()