from numpy import median
f = open("./reportFile.csv").readlines()


def findMZinmzDic(mz, charge, mzdic):
    if (mz, charge) in mzdic:
        return mzdic[mz, charge]
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

print(repDic)

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

combDic = {}
for group in newDic:
    mz, charge = newDic[group][:2]
    if not findMZinmzDic(mz, charge, combDic):
        combDic[mz, charge] = [newDic[group]]
    else:
        matched = findMZinmzDic(mz, charge, combDic)
        combDic[matched].append(newDic[group])

print(combDic)
for wd in combDic:
    combDic[wd] = sorted(combDic[wd], key= lambda x: x[2], reverse=True)

print(combDic)
for wd in combDic:
    if len(combDic[wd]) <= 2:
        pass
    else:
        combDic[wd] = combDic[wd][:2]
print(len(combDic))

scanTargetList = []
for wd in combDic:
    for info in combDic[wd]:
        bg, ed = info[3]
        scanTargetList.extend(list(range(bg, ed+1)))

scanTargetList.sort()
print(scanTargetList)


