import os
from detect_isotopic import detectIsotopic

# os.chdir(r'G:\20190304_Lumos\PEPTIDES')

def generate_mass_range(num, delta_ppm):
    delta = num * delta_ppm / 1000000
    return num - delta, num + delta


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


def readMS1file(openedMS1, scanNum, scanRange):
    eluteMS1dic = {}
    lowScanNum = scanNum - scanRange
    upScanNum = scanNum + scanRange
    print(lowScanNum, upScanNum)
    i = 0
    while i < len(openedMS1):
        if openedMS1[i][0] != "S":
            i += 1
        else:
            currentScan = int(openedMS1[i].split("\t")[1])
            if currentScan < lowScanNum:
                i += 1
            else:
                if currentScan > upScanNum:
                    break
                else:
                    currentRT = openedMS1[i+2].strip().split("\t")[-1]
                    currentMzIntensDic = {}
                    p = i + 5
                    while p < len(openedMS1):
                        if openedMS1[p][0] == 'S':
                            break
                        else:
                            [mz, intens] = openedMS1[p].strip().split(" ")
                            currentMzIntensDic[float(mz)] = float(intens)
                        p += 1
                    eluteMS1dic[currentScan] = [currentRT, currentMzIntensDic]
                    i = p

    return eluteMS1dic


def calcEluteArea(eluteMS1dic, mz):
    wantedMZeluteDic = {}
    for scan in eluteMS1dic:
        currRT, currMzInDic = eluteMS1dic[scan]
        chargedMZdic = detectIsotopic(currMzInDic)
        currMzList = sorted(list(chargedMZdic.keys()))
        findBool, matchedMZ = findOneMZ1(mz, currMzList, 10)
        if findBool:
            wantedMZeluteDic[scan] = [currRT, currMzInDic[matchedMZ]]
        else:
            wantedMZeluteDic[scan] = [currRT, 0]
    return wantedMZeluteDic


def main():
    f = open('VR_6_PHGO.ms1', 'r').readlines()
    eluteRangeDic = readMS1file(f, 4930, 100)
    wantedMZeluteDic = calcEluteArea(eluteRangeDic, 387.2305)
    print(wantedMZeluteDic)


if __name__ == "__main__":
    main()

#scan1dic = readMS1file(f, 2, 0)[2][1]
# print(scan1dic)
#scan1MZlist = sorted(list(scan1dic.keys()))
#for key in scan1MZlist:
    #print(key, scan1dic[key])