# coding = utf-8
import os
from numpy import median

#os.chdir(r"F:\MS_DATA_STORAGE\20190819\multiNCE\incul_random")

def findMZinmzDic(mz, charge, mzdic):
    if (mz, charge) in mzdic:
        return mz, charge
    else:
        mzList = sorted([x[0] for x in list(mzdic.keys()) if x[1] == charge])
        if mzList == []:
            return False
        else:
            deltaMass = mz * 10 / 1000000
            lowMZ, upMZ = mz - deltaMass, mz + deltaMass
            if lowMZ > mzList[-1]:
                return False
            else:
                i = 0
                while i < len(mzList):
                    if upMZ < mzList[i]:
                        return False
                    elif lowMZ <= mzList[i]:
                        if (mzList[i], charge) in mzdic:
                            return mzList[i], charge
                        else:
                            return False
                    else:
                        i += 1


def groupScans(openedFL):
    f = openedFL
    repDic = {}
    i = 1
    while i < len(f):
        print(i)
        scan = int(f[i].strip().split(",")[0])
        endGroup = scan + 9
        p = i
        infoList = []
        while p < i +10:
            infoList.append(f[p])
            p += 1
        infoList = sorted(infoList, key = lambda x: int(x.split(",")[3]))
        repDic[scan, endGroup] = infoList
        i = p 
    return repDic


def groupMZchargeCenter(grouped_repDic):
    repDic = grouped_repDic
    newDic = {}
    for group in repDic:
        infoList = repDic[group]
        commBYratioList = []
        reportNum = 0
        mzList = []
        chargeList = []
        for info in infoList:
            wdlist = info.strip().split(",")
            if wdlist[1] not in mzList:
                mzList.append(wdlist[1])
            if wdlist[2] not in chargeList:
                chargeList.append(wdlist[2])
            commBYratioList.append(float(wdlist[6]))
            reportNum += int(wdlist[9])
        if len(mzList) != 1 or len(chargeList) != 1:
            print("wrong" + group)
        else:
            mz = float(mzList[0])
            charge = chargeList[0]
            commRatioMedia =  median(commBYratioList)
            reportNumRatio = reportNum/20
            score = commRatioMedia*0.7 + reportNumRatio * 0.3
            if commRatioMedia > 0.5 and reportNumRatio > 0.4:
                newDic[group] = [mz, charge, score, group, commRatioMedia, reportNumRatio]
    return newDic


def removeReplicate(newDic):
    combDic = {}
    for group in newDic:
        mz, charge = newDic[group][:2]
        if not findMZinmzDic(mz, charge, combDic):
            combDic[mz, charge] = [newDic[group]]
        else:
            matchedmz, charge = findMZinmzDic(mz, charge, combDic)
            combDic[matchedmz, charge].append(newDic[group])

    for wd in combDic:
        combDic[wd] = sorted(combDic[wd], key= lambda x: x[2], reverse=True)

    for wd in combDic:
        if len(combDic[wd]) <= 2:
            pass
        else:
            combDic[wd] = combDic[wd][:2]
    return combDic


def filterResult(combDic, openedFL):
    f = openedFL
    scanTargetList = []
    for wd in combDic:
        for info in combDic[wd]:
            bg, ed = info[3]
            scanTargetList.extend(list(range(bg, ed+1)))

    scanTargetList.sort()

    b = open("resultFilter.csv", 'w')
    b.write(f[0])
    for line in f[1:]:
        if int(line.strip().split(",")[0]) in scanTargetList:
            b.write(line)
    b.close()


def loadtypeColDic():
    typeColumnDic = {14: "commBYbeta", 13:"commBYalpha", 12: "interPXratio", \
        11:"rep_ints_alpha", 10: "rep_ints_beta", 9: "pair_num", \
        7: "XlinkBYratio", 6: "commBYratio"}
    return typeColumnDic


def extractInfo(column, filename, openedFilterFile):
    ft = openedFilterFile
    c = open(filename + ".csv", 'w')
    combyNCEdic = {}
    for line in ft[1:]:
        lineList = line.strip().split(",")
        tagt = lineList[column]
        nce = lineList[3]
        if nce not in combyNCEdic:
            combyNCEdic[nce] = [tagt]
        else:
            combyNCEdic[nce].append(tagt)

    nceList = sorted(list(combyNCEdic.keys()))
    print(nceList)
    for i in range(len(combyNCEdic["22"])):
        wlist = []
        for nce in nceList:
            wlist.append(combyNCEdic[nce][i])
        c.write(",".join(wlist)+"\n")

    c.close()


def main():
    typeColumnDic = loadtypeColDic()
    f = open("./reportFile_commBY.csv").readlines()
    groupScanDic = groupScans(f)
    mzchgCentDic = groupMZchargeCenter(groupScanDic)
    finalDic = removeReplicate(mzchgCentDic)
    filterResult(finalDic, f)
    ft = open("resultFilter.csv", 'r').readlines()
    for colmun in typeColumnDic:
        filename = typeColumnDic[colmun]
        extractInfo(colmun, filename, ft)
    

if __name__ == "__main__":
    main()