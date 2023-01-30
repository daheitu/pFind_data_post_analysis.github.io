import os
import random
from detect_isotopic import detectIsotopic


# os.chdir(r'G:\6peptide_PhGO_HF_20190223')

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
    if k == len(ms1MzList):
        return False, -1


def detectPepFeat(openedMS1, intensCutoff):
    finalPepFeatDic = {}
    i = 0
    while i < len(openedMS1):
        if openedMS1[i][0] != "S":
            i += 1
        else:
            currentScan = int(openedMS1[i].split("\t")[1])
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
            currChgedIsoDic = detectIsotopic(currentMzIntensDic)
            for mz in currChgedIsoDic:
                chrg = currChgedIsoDic[mz][0]
                intens = currentMzIntensDic[mz]
                if chrg == 1:
                    continue
                else:
                    if intens < intensCutoff:
                        continue
                    else:
                        if len(finalPepFeatDic):
                            currFanalList = sorted(list(finalPepFeatDic.keys()))
                            [findBool, matchedMZ] = findOneMZ1(mz, currFanalList, 10)
                            if findBool == False:
                                finalPepFeatDic[mz] = [mz, chrg, intens, currentScan]
                            else:
                                if intens < finalPepFeatDic[matchedMZ][2]:
                                     pass
                                else:
                                    del finalPepFeatDic[matchedMZ]
                                    finalPepFeatDic[mz] = [mz, chrg, intens, currentScan]
                        else:
                            finalPepFeatDic[mz] = [mz, chrg, intens, currentScan]

            i = p
    return finalPepFeatDic


def readMS1file(openedMS1, scanNum, scanRange):
    eluteMS1dic = {}
    lowScanNum = scanNum - scanRange
    upScanNum = scanNum + scanRange
    # print(lowScanNum, upScanNum)
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


def detecLocalBaseLine(eluteMS1dic, mz):
    baseIntensList = []
    for scan in eluteMS1dic:
        currRT, currMzInDic = eluteMS1dic[scan]
        currMzList = sorted(list(currMzInDic.keys()))
        findBool, matchedMZ = findOneMZ1(mz, currMzList, 10)
        if findBool:
            baseIntensList.append(currMzInDic[matchedMZ])
    baseIntens = sum(baseIntensList)/len(baseIntensList)
    return baseIntens


def detectGlobalBaseLine(openedMS1, mz, randNum):
    for line in openedMS1[::-1]:
        if line[0] == "S":
            totalScan = int(line.split("\t")[1])
            break
        else:
            continue
    
    intensList = []
    for i in range(0, randNum):
        seed = random.randint(1, totalScan)
        eluteRangeDic = readMS1file(openedMS1, seed, 100)
        intensCutoff = detecLocalBaseLine(eluteRangeDic, mz)
        intensList.append(intensCutoff)
    
    return round(sum(intensList)/len(intensList), 1)


def calcEluteArea(eluteMS1dic, mz, intes, peakbase):
    intensCutoff = intes * peakbase
    wantedMZeluteDic = {}
    for scan in eluteMS1dic:
        currRT, currMzInDic = eluteMS1dic[scan]
        chargedMZdic = detectIsotopic(currMzInDic)
        currMzList = sorted(list(chargedMZdic.keys()))
        findBool, matchedMZ = findOneMZ1(mz, currMzList, 20)
        if findBool:
            wantedMZeluteDic[scan] = [currRT, currMzInDic[matchedMZ]]
        else:
            wantedMZeluteDic[scan] = [currRT, 0]
    
    scanList = sorted(list(wantedMZeluteDic.keys()))
    areaTotal = 0
    for i in range(len(scanList)-1):
        currRT, CurrIntens = wantedMZeluteDic[scanList[i]]
        nextRT, nextIntens = wantedMZeluteDic[scanList[i + 1]]
        if max(CurrIntens, nextIntens) < intensCutoff:
            pass
        else:
            area = (CurrIntens + nextIntens) * (float(nextRT) -float(currRT)) /2
        
            areaTotal += round(area, 2)
    return areaTotal


def main():
    peakbase = 0.1
    for fl in os.listdir(os.getcwd()):
        if fl[-4:] == ".ms1":
            print("the current ms file is " + fl)
            f = open(fl, 'r').readlines()
            repName = fl[:-4] + "pep_list_perc10.txt"
            b = open(repName, 'w')
            b.write("\t".join(["mz", "charge", "Intensity", "scanNum", "Area"]) + "\n")
            baseIntens = detectGlobalBaseLine(f, 462.14658, 5)
            finalPepFeatDic = detectPepFeat(f, baseIntens)
            for mz in finalPepFeatDic:
                scan = finalPepFeatDic[mz][-1]
                intes = finalPepFeatDic[mz][-2]
                eluteRangeDic = readMS1file(f, scan, 100)
                wantedMZeluteArea = calcEluteArea(eluteRangeDic, mz, intes, peakbase)
                finalPepFeatDic[mz].append(wantedMZeluteArea)
                writeLine = [str(ele) for ele in finalPepFeatDic[mz]]
                b.write("\t".join(writeLine) + "\n")
            
            b.close()


if __name__ == "__main__":
    main()
