import os
from copy import deepcopy
import numpy as np


# 读取ms3文件，获得谱图号以及其对应的mz，charge，队列能量等信息
def readms3Info(ms3Path):
    ms3InfoDic = {}
    f = open(ms3Path, 'r').readlines()
    i = 0
    while i < len(f):
        if f[i][0] != "S":
            i += 1
        else:
            scanNum = int(f[i].strip().split("\t")[1])
            precMS2mz, nceMS2 = f[i+6].strip().split(" ")[7].split("@")
            precMS3mz, nceMS3 = f[i+6].strip().split(" ")[8].split("@")
            precMS2mz = float(precMS2mz)
            precMS3mz = float(precMS3mz)
            actTypeMS2= "".join([ele for ele in nceMS2 if ele.isalpha()])
            engMS2 = int(float(nceMS2[nceMS2.find(actTypeMS2)+len(actTypeMS2):]))
            actTypeMS3= "".join([ele for ele in nceMS3 if ele.isalpha()])
            engMS3 = int(float(nceMS3[nceMS3.find(actTypeMS3)+len(actTypeMS3):]))
            scanPrecMS2 = int(f[i+7].strip().split("\t")[2])
            charge = int(f[i+9].strip().split("\t")[1])
            ms3InfoDic[scanNum] = [scanNum, precMS3mz, charge, actTypeMS3, engMS3, scanPrecMS2, precMS2mz, actTypeMS2, engMS2]
            i += 10
    return ms3InfoDic


# 根据ms3的信息，将来自同一ms2 的归并到ms2——triggered字典里
def getLinageInfo(ms3InfoDic):
    ms2TriggerMS3dic = {}
    for scan in ms3InfoDic:
        scanPrecMS2 = ms3InfoDic[scan][5]
        precMS3mz = ms3InfoDic[scan][1]
        precMS3charge = ms3InfoDic[scan][2]
        if scanPrecMS2 not in ms2TriggerMS3dic:
            ms2TriggerMS3dic[scanPrecMS2] = [(precMS3mz, precMS3charge, scan)]
        else:
            ms2TriggerMS3dic[scanPrecMS2].append((precMS3mz, precMS3charge, scan))
    
    return ms2TriggerMS3dic


# 遍历ms2文件，返回总的ms2的数目以及每个被触发三级的ms2的母离子质量和价态
def readms2Info(ms2path, ms2TriggerMS3dic):
    f = open(ms2path, 'r')
    ms2f = f.readlines()
    totalNumOfms2 = 0
    i = 0
    while i < len(ms2f):
        if ms2f[i][0] != "S":
            i += 1
        else:
            totalNumOfms2 += 1
            scan = int(ms2f[i].strip().split("\t")[1])
            mz = float(ms2f[i].strip().split("\t")[-1])
            if scan not in ms2TriggerMS3dic:
                i += 10
            else:
                charge = int(ms2f[i+9].strip().split('\t')[1])
                massH1 = float(ms2f[i+9].strip().split('\t')[-1])
                ms2TriggerMS3dic[scan].insert(0, (scan, mz, charge, massH1))
                i += 10
    return totalNumOfms2, ms2TriggerMS3dic



def dealstrtuple(pairedInfo):
    pairedInfoList = []
    if ";" not in pairedInfo:
        print("wrong")
    else:
        for key in pairedInfo.split(";"):
            [scan1, scan2] = key[1:-1].split(",")
            scan1 = int(scan1); scan2 = int(scan2)
            pairedInfoList.append((scan1, scan2))
    return pairedInfoList


# 把ms2 triggered 信息写入文件
def writeMs2ToMs3Info(ms2TriggerMS3dic, repName):
    b = open(repName, 'w')
    for ms2scan in ms2TriggerMS3dic:
        writeList = [str(ms2scan)]
        for ion in ms2TriggerMS3dic[ms2scan]:
            writeList.append(str(ion))
        b.write("\t".join(writeList) + "\n")
    b.close()



def compareTwoNum(therMass, realMass, tol_ppm):
    deltaPPM = abs(realMass - therMass)/therMass*1000000
    if deltaPPM < tol_ppm :
        return True, round(deltaPPM, 2)
    else:
        return False, round(deltaPPM, 2)


def detectChargeFromMZ(mz1, mz2):
    deltaMzChrageDic = {}
    for i in range(1,6):
        deltaMZ = 31.9720/i
        deltaMzChrageDic[deltaMZ] = i
    
    deltaMzList = deltaMzChrageDic.keys()
    deltaMZ1MZ2_colse = 0
    deltaMZ1MZ2 = abs(mz1-mz2)
    if deltaMZ1MZ2 > 31.9720:
        return 0
    else:
        for mz in deltaMzList:
            if abs(mz - deltaMZ1MZ2) < abs(deltaMZ1MZ2_colse-deltaMZ1MZ2):
                deltaMZ1MZ2_colse = mz
        if deltaMZ1MZ2_colse == 0:
            chg = 0
        else:
            chg = deltaMzChrageDic[deltaMZ1MZ2_colse]
        return chg


def detectCharge(mz1, charge1, mz2, charge2):
    if charge1 == charge2:
        if charge1 == 0:
            print("both charge is 0")
            chg = detectChargeFromMZ(mz1, mz2)
        else:
            chg = charge1
    else:
        if charge1 == 0:
            chg = charge2
        elif charge2 == 0:
            chg = charge1
        else:
            chg = 0
    return chg


def judgePair(feat1Info, feat2Info):
    mz1, charge1 = feat1Info[:-1]
    mz2, charge2 = feat2Info[:-1]
    chg = detectCharge(mz1, charge1, mz2, charge2)
    
    if chg == 0:
        return False, abs(mz1-mz2)
    else:
        matchBool, delta = compareTwoNum(min(mz1, mz2) + 31.972/chg, max(mz1,mz2), 30)
        if matchBool:
            return True, str(delta) + "ppm" 
        else:
            return False, str(delta) + "ppm" 


def findPairfromList(featList):
    pairedFeatList = []
    unpairedFeatList = []
    lenFeatList = len(featList)
    if lenFeatList < 2 or lenFeatList > 4:
        print("featList format is wrong")
    else:
        for i in range(lenFeatList):
            for j in range(i+1, lenFeatList):
                if judgePair(featList[i], featList[j])[0]:
                    pairedFeatList.append((featList[i], featList[j]))
                    break
                else:
                    continue
            if j == lenFeatList -1:
                unpairedFeatList.append(featList[i])
    return pairedFeatList, unpairedFeatList
    
                 
def linkJudgetwoMZ(preMHpair1, preMHpair2, precuMH):
    sumAB = preMHpair1+preMHpair2
    longShpairBool = compareTwoNum(precuMH+1.0078-18.0105, sumAB, 10)[0]
    if longShpairBool:
        return True
    else:
        shortShpairBool = compareTwoNum(precuMH+1.0078-18.0105-31.97207, sumAB, 10)[0]
        if shortShpairBool:
            return True
        else:
            longLpairBool = compareTwoNum(precuMH+1.0078-18.0105+31.97207, sumAB, 10)[0]
            if longLpairBool:
                return True
            else:
                return False
        

def monoJudgeOneMZ(ms3Mass, precMass):
    islongBool = compareTwoNum(precMass-176.015+85.983, ms3Mass, 10)[0]
    if islongBool:
        return True
    else:
        isshortBool = compareTwoNum(precMass-176.015+54.011, ms3Mass, 10)[0]
        if isshortBool:
            return True
        else:
            return False


def staticMs2TriggerMS3dic(ms2TriggerMS3dic):
    finalDic = {}
    unpairedMS2 = 0
    unpairedlenDic = {}
    for scan in ms2TriggerMS3dic:
        precuMH =  ms2TriggerMS3dic[scan][0][-1]
        numOfms3forMS2 = len(ms2TriggerMS3dic[scan]) - 2
        if numOfms3forMS2 % 2 == 1:
            unpairedMS2 += 1
            if numOfms3forMS2 not in unpairedlenDic:
                unpairedlenDic[numOfms3forMS2] = 1
            else:
                unpairedlenDic[numOfms3forMS2] += 1
            unpairedMS2 += 1
        
        if numOfms3forMS2 == 1:
            ismonoBool = monoJudgeOneMZ(ms2TriggerMS3dic[scan][-2], precuMH)
            islinkedBool = False
        elif numOfms3forMS2 == 2:
            feat1Info, feat2Info = ms2TriggerMS3dic[scan][1:-1]
            if judgePair(feat1Info, feat2Info):
                ismonoBool = monoJudgeOneMZ(feat1Info, precuMH)
                islinkedBool = False
            else:
                print("scan " + str(scan) + " not paired in two")
                ismonoBool = monoJudgeOneMZ(feat1Info, precuMH) or monoJudgeOneMZ(feat2Info, precuMH)
                islinkedBool = False
        elif numOfms3forMS2 == 3:
            pairedfeaList, unpairedfeaList = findPairfromList(ms2TriggerMS3dic[scan][1:-1])
            if len(pairedfeaList) == 0:
                print("scan " + str(scan) + " not paired in three")
            else:
                pair1Info = pairedfeaList[0][0]
                ismonoBool = monoJudgeOneMZ(pair1Info, precuMH)
                if not ismonoBool:
                    islinkedBool = linkJudgetwoMZ(pair1Info, unpairedfeaList[0], precuMH)
                else:
                    islinkedBool = False
        else:
            pairedfeaList, unpairedfeaList = findPairfromList(ms2TriggerMS3dic[scan][1:])
            if len(pairedfeaList) in [0, 1]:
                print("scan " + str(scan) + " not paired in four")
            else:
                ismonoBool = False
                for pair in pairedfeaList:
                    pair11Info = pair[0]
                    if monoJudgeOneMZ(pair11Info, precuMH):
                        ismonoBool = True
                if ismonoBool:
                    islinkedBool = "not sure"
                else:
                    pair1Info = pairedfeaList[0][0]
                    pair2Info = pairedfeaList[1][0]
                    islinkedBool = linkJudgetwoMZ(pair1Info, pair2Info, precuMH)
        repoList = ms2TriggerMS3dic[scan]
        repoList.insert(1, ismonoBool)
        repoList.insert(2, islinkedBool)
        finalDic[scan] = repoList
    unpairedRatio = unpairedMS2/len(ms2TriggerMS3dic)
    return unpairedRatio, unpairedlenDic, finalDic


def staticUnmatchMs3(ms2TriggerMS3dic, repName):
    numOfUnpairedMS3 = 0
    for scan in ms2TriggerMS3dic:
        numOfTriggedMS3 = len(ms2TriggerMS3dic[scan])-2
        pairedInfo = ms2TriggerMS3dic[scan][-1]
        if numOfTriggedMS3 == 1:
            numOfUnpairedMS3 += 1
        elif numOfTriggedMS3 == 2:
            if pairedInfo == "None":
                numOfUnpairedMS3 += 2
        elif numOfTriggedMS3 == 3:
            if pairedInfo == "None":
                numOfUnpairedMS3 +=3
            else:
                numOfUnpairedMS3 += 1
        elif numOfTriggedMS3 == 4:
            if ";" not in pairedInfo:
                numOfUnpairedMS3 += 2
    return numOfUnpairedMS3


def calPairedIonDeltaPPM(ms2TriggerMS3dic):
    repDic = {}
    for ms2scan in ms2TriggerMS3dic:
        numOfTriggedMS3 = len(ms2TriggerMS3dic[ms2scan])-2
        pairedInfo = ms2TriggerMS3dic[ms2scan][-1]
        featIonInfoList = ms2TriggerMS3dic[ms2scan][1:-1]
        pairDeltaDic = {}
        if pairedInfo != "None":
            if ';' not in pairedInfo:
                pairedFeatIonList = []
                for scan in pairedInfo:
                    for featIon in featIonInfoList:
                        if scan in featIon:
                            pairedFeatIonList.append(featIon)
                feat1Info, feat2Info = pairedFeatIonList
                deltappm = judgePair(feat1Info, feat2Info)[1]
                chargePair = (feat1Info[1], feat2Info[1])
                pairDeltaDic[pairedInfo] = [deltappm, chargePair]

            else:
                pairedInfoList = dealstrtuple(pairedInfo)
                #print(pairedInfoList)
                for pairedInfo in pairedInfoList:
                    pairedFeatIonList = []
                    for scan in pairedInfo:
                        for featIon in featIonInfoList:
                            if scan in featIon:
                                pairedFeatIonList.append(featIon)
                    feat1Info, feat2Info = pairedFeatIonList
                    deltappm = judgePair(feat1Info, feat2Info)[1]
                    chargePair = (feat1Info[1], feat2Info[1])
                    pairDeltaDic[pairedInfo] = [deltappm, chargePair]
        repDic[ms2scan] = pairDeltaDic
    
    return repDic


def calunmatcheddelta(ms2TriggerMS3dic, fltoWT):
    b = open(fltoWT, "w")
    repDic = {}
    numWrongCal = 0
    for scan in ms2TriggerMS3dic:
        repList = deepcopy(ms2TriggerMS3dic[scan])
        featInfoList = ms2TriggerMS3dic[scan][1:]
        if len(featInfoList) == 1:
            ms2TriggerMS3dic[scan].append("None")
        elif len(featInfoList) == 2:
            feat1Info, feat2Info = featInfoList
            mactchBool, delta = judgePair(feat1Info, feat2Info)
            if not mactchBool:
                numWrongCal += 2
                repList.append(delta)
                repList.insert(0, scan)
                b.write("\t".join([str(ele) for ele in repList])+"\n")
                ms2TriggerMS3dic[scan].append("None")
            else:
                pairedInfo = (feat1Info[-1], feat2Info[-1])
                ms2TriggerMS3dic[scan].append(pairedInfo)
        elif len(featInfoList) == 4:
            feat1Info, feat2Info = featInfoList[:2]
            feat3Info, feat4Info = featInfoList[2:]
            pairedInfoList = []
            isappendBool = False
            for (pair1, pair2) in [(feat1Info, feat2Info), (feat3Info, feat4Info)]:
                mactchBool, delta = judgePair(pair1, pair2)
                if not mactchBool:
                    numWrongCal += 2
                    repList.insert(0, scan)
                    repList.append((pair1[-1], pair2[-1]))
                    repList.append(delta)
                    isappendBool = True     
                else:
                    pairedInfo = (pair1[-1], pair2[-1])
                    pairedInfoList.append(pairedInfo)
            if len(pairedInfoList) == 0:
                ms2TriggerMS3dic[scan].append("None")
            elif len(pairedInfoList) == 1:
                ms2TriggerMS3dic[scan].append(pairedInfo)
            else:
                ms2TriggerMS3dic[scan].append(";".join([str(ele) for ele in pairedInfoList]))

            if isappendBool:
                b.write("\t".join([str(ele) for ele in repList])+"\n")
        elif len(featInfoList) == 3:
            feat1Info, feat2Info, feat3Info = featInfoList[:3]
            mactchBool12, delta12 = judgePair(feat1Info, feat2Info)
            mactchBool23, delta23 = judgePair(feat3Info, feat2Info)
            if not mactchBool12 and not mactchBool23:
                numWrongCal += 3
                repList.insert(0, scan)
                repList.append((feat1Info[-1], feat2Info[-1]))
                repList.append(delta12)
                repList.append((feat2Info[-1], feat3Info[-1]))
                repList.append(delta23)
                b.write("\t".join([str(ele) for ele in repList])+"\n")
                ms2TriggerMS3dic[scan].append("None")
                    #ms2TriggerMS3dic[scan].append(pairedInfo)
            elif mactchBool12:
                pairedInfo = (feat1Info[-1], feat2Info[-1])
                ms2TriggerMS3dic[scan].append(pairedInfo)
            elif mactchBool23:
                pairedInfo = (feat2Info[-1], feat3Info[-1])
                ms2TriggerMS3dic[scan].append(pairedInfo)
    b.close()
    return numWrongCal


def pickupChargeZeroScan(ms2TriggerMS3dic, repFl):
    numOFnondetect = 0
    b = open(repFl, 'w')
    for scan in ms2TriggerMS3dic:
        featInfoList = ms2TriggerMS3dic[scan][1:-1]
        pairedInfo = ms2TriggerMS3dic[scan][-1]
        contain0chrgeBool = False
        for featInfo in featInfoList:
            charge = featInfo[1]
            scan = featInfo[2]
            if charge == 0:
                contain0chrgeBool = True
                if str(scan) in pairedInfo:
                    numOFnondetect += 2
                else:
                    numOFnondetect += 1
                break
            else:
                continue
        if contain0chrgeBool:
            writeList = [str(scan)]
            for featInfo in featInfoList:
                writeList.append(str(featInfo))
            b.write("\t".join(writeList) + "\n")
    b.close()
    return numOFnondetect


def writeDeltaPPM(paireDELTADic, repName):
    b = open(repName, 'w')
    chargeZeroDeltaList = []
    otherDelatList = []
    for scan in paireDELTADic:
        pairDeltaInnerDic = paireDELTADic[scan]
        for pair in pairDeltaInnerDic:
            chargePair = pairDeltaInnerDic[pair][1]
            dletappm = pairDeltaInnerDic[pair][0]
            wlist = [str(scan), str(pair), str(chargePair), str(dletappm)]
            b.write("\t".join(wlist)+'\n')
            if 0 in chargePair:
                chargeZeroDeltaList.append(dletappm)
            else:
                otherDelatList.append(dletappm)
    b.close()
    return chargeZeroDeltaList, otherDelatList


def findScanInfo(scan, featIonInfoList):
    for featIon in featIonInfoList:
        if featIon[2] == scan:
            return featIon
        else:
            continue


def detectChargeZero(ms2TriggerMS3dic):
    unpairedChargeZeroDic = {}
    pairedChargeZeroDic = {}
    for scan in ms2TriggerMS3dic:
        numOfTriggedMS3 = len(ms2TriggerMS3dic[scan])-2
        pairedInfo = ms2TriggerMS3dic[scan][-1]
        featIonInfoList = ms2TriggerMS3dic[scan][1:-1]
        scan_info_dic = {}
        if ";" in pairedInfo:
            pairedInfoList = dealstrtuple(pairedInfo)
        else:
            pairedInfoList = [pairedInfo]
        if pairedInfo != "None":
            continue
        else:
            for featIon in featIonInfoList:
                if featIon[1] == 0:
                    for pairInfo in pairedInfoList:
                        if featIon[2] in pairInfo:
                            tgtFeatInfoList = []
                            for scan in pairInfo:
                                tgtfeatIon = findScanInfo(scan, featIonInfoList)
                                tgtFeatInfoList.append(tgtfeatIon)
                            scan_info_dic[pairInfo] = tgtFeatInfoList
                        else:
                            scan_info_dic[scan] = findScanInfo(scan, featIonInfoList)
                else:
                    continue
                    #if str(featIon[2]) not in str(pairedInfo):

            if ';' not in pairedInfo:
                pairedFeatIonList = []
                for scan in pairedInfo:
                    for featIon in featIonInfoList:
                        if scan in featIon:
                            pairedFeatIonList.append(featIon)
                feat1Info, feat2Info = pairedFeatIonList
                deltappm = judgePair(feat1Info, feat2Info)[1]
                pairDeltaDic[pairedInfo] = deltappm
            else:
                pairedInfoList = dealstrtuple(pairedInfo)
                #print(pairedInfoList)
                for pairedInfo in pairedInfoList:
                    pairedFeatIonList = []
                    for scan in pairedInfo:
                        for featIon in featIonInfoList:
                            if scan in featIon:
                                pairedFeatIonList.append(featIon)
                    feat1Info, feat2Info = pairedFeatIonList
                    deltappm = judgePair(feat1Info, feat2Info)[1]
                    pairDeltaDic[pairedInfo] = deltappm
        repDic[scan] = pairDeltaDic


def writedeltaList(wlist, wflname):
    b = open(wflname, 'w')
    for ele in wlist:
        b.write(ele[:-3]+"\n")
    b.close()


def stacMS2MS3(ms2TriggerMS3dic):
    repDic = {}
    for ms2scan in ms2TriggerMS3dic:
        pairState = ms2TriggerMS3dic[ms2scan][-1]
        if pairState == "None":
            state = "None"
        elif ";" not in pairState:
            state = "onePair"
        else:
            state = "twoPair"
        numOFMS3 = len(ms2TriggerMS3dic[ms2scan]) -2
        if state not in repDic:
            repDic[state] = {}
            if  numOFMS3 not in repDic[state]:
                repDic[state][numOFMS3] = 1
            else:
                repDic[state][numOFMS3] += 1
        else:
            if  numOFMS3 not in repDic[state]:
                repDic[state][numOFMS3] = 1
            else:
                repDic[state][numOFMS3] += 1
        
    return repDic
    

def calIonMH(mz, charge):
    mz = float(mz)
    chr = int(charge)
    ionMH = (mz-1.0078)*chr + 1.0078
    
    return round(ionMH, 5)


def findPairInfo(ms3pair, ms3infolist):
    for ms3info in ms3infolist:
        if 0 in ms3info:
            continue
        else:
            if int(ms3info[-1]) in ms3pair:
                mzMs3, chrms3 = ms3info[:-1]
                return calIonMH(mzMs3, chrms3)
            else:
                continue
                

def stac2pairinfo(ms2TriggerMS3dic):
    num_xlink_dic = {}
    num_mono_dic = {}
    others_dic = {}
    for ms2can in ms2TriggerMS3dic:
        pairState = ms2TriggerMS3dic[ms2can][-1]
        if ";" not in pairState:
            continue
        else:
            precuMH = ms2TriggerMS3dic[ms2can][0][-1]
            pair1, pair2 = dealstrtuple(pairState)
            ms3infolist = ms2TriggerMS3dic[ms2can][1:-1]
            preMHpair1 = findPairInfo(pair1, ms3infolist)
            preMHpair2 = findPairInfo(pair2, ms3infolist)
            if linkJudgetwoMZ(preMHpair1, preMHpair2, precuMH):
                num_xlink_dic[ms2can] = ms2TriggerMS3dic[ms2can]
            else:
                boolcontMono = False
                for msofms3 in [preMHpair1, preMHpair2]:
                    if monoJudgeOneMZ(msofms3, precuMH):
                        boolcontMono = True
                        break
                    else:
                        continue
                if boolcontMono:
                    num_mono_dic[ms2can] = ms2TriggerMS3dic[ms2can]
                else:
                    others_dic[ms2can] = ms2TriggerMS3dic[ms2can]
    return num_xlink_dic, num_mono_dic, others_dic


def stacOnePairinfo(ms2TriggerMS3dic):
    monoDic = {}
    otherdic = {}
    xlinkDic = {}
    for ms2scan in ms2TriggerMS3dic:
        pairState = ms2TriggerMS3dic[ms2scan][-1]
        if type(pairState) != tuple:
            continue
        else:
            precuMH = ms2TriggerMS3dic[ms2scan][0][-1]
            ms3InfoList = ms2TriggerMS3dic[ms2scan][1:-1]
            if len(ms3InfoList) == 2:
                preMHms3 = findPairInfo(pairState, ms3InfoList)
                if monoJudgeOneMZ(preMHms3, precuMH):
                    monoDic[ms2scan] = ms2TriggerMS3dic[ms2scan]
                else:
                    otherdic[ms2scan] = ms2TriggerMS3dic[ms2scan]
            elif len(ms3InfoList) == 3:
                preMHms3 = findPairInfo(pairState, ms3InfoList)
                for ms3info in ms3InfoList:
                    scan = ms3info[-1]
                    if scan not in pairState:
                        if ms3info[1] == 0:
                            if monoJudgeOneMZ(preMHms3, precuMH):
                                monoDic[ms2scan] = ms2TriggerMS3dic[ms2scan]
                            else:
                                otherdic[ms2scan] = ms2TriggerMS3dic[ms2scan]
                        else:
                            aloneMH = calIonMH(ms3info[0], ms3info[1])
                            if linkJudgetwoMZ(aloneMH, preMHms3, precuMH):
                                xlinkDic[ms2scan] = ms2TriggerMS3dic[ms2scan]
                            else:
                                boolcontMono = False
                                for msofms3 in [aloneMH, preMHms3]:
                                    if monoJudgeOneMZ(msofms3, precuMH):
                                        boolcontMono = True
                                        break
                                    else:
                                        continue
                                if boolcontMono:
                                    monoDic[ms2scan] = ms2TriggerMS3dic[ms2scan]
                                else:
                                    otherdic[ms2scan] = ms2TriggerMS3dic[ms2scan]
    return xlinkDic, monoDic, otherdic


def main():
    ms3InfoDic = readms3Info("./DSSO_CID_MS2_CID_MS3_T1.ms3")
    numTotalMS3 = len(ms3InfoDic)
    ms2TriggerMS3dic = getLinageInfo(ms3InfoDic)
    totalms2, ms2TriggerMS3dic = readms2Info("./DSSO_CID_MS2_CID_MS3_T1.ms2", ms2TriggerMS3dic)
    writeMs2ToMs3Info(ms2TriggerMS3dic, "ms2_ms3_linage.txt")
    # triggedNumMs2 = len(ms2TriggerMS3dic)
    # triggedRatio = triggedNumMs2/totalms2
    
    # numWrongCal = calunmatcheddelta(ms2TriggerMS3dic, "unmatchedScans.txt")
    # ratioWrongCal = numWrongCal/numTotalMS3
    # writeMs2ToMs3Info(ms2TriggerMS3dic, "paired_info.txt")
    # numOFnondetect = pickupChargeZeroScan(ms2TriggerMS3dic, "charge_zero.txt")
    # ratiowrongcharge = numOFnondetect/numTotalMS3
    # print(ratioWrongCal, ratiowrongcharge)
    # paireDELTADic = calPairedIonDeltaPPM(ms2TriggerMS3dic)
    # zeroMed, otherMed = writeDeltaPPM(paireDELTADic, "deltappm.txt")
    # #print(zeroMed, otherMed)
    # writedeltaList(zeroMed, 'zeroMed.txt')
    # writedeltaList(otherMed, 'otherMed.txt')
    # print(stacMS2MS3(ms2TriggerMS3dic))
    # xlinkInTwoPDic, monoInTwoPDic, otherInTwoPDic = stac2pairinfo(ms2TriggerMS3dic)
    # print(len(xlinkInTwoPDic), len(monoInTwoPDic), len(otherInTwoPDic))
    # xlink2deltaDic = calPairedIonDeltaPPM(xlinkInTwoPDic)
    # mono2deltaDic = calPairedIonDeltaPPM(monoInTwoPDic)
    # other2deltaDic = calPairedIonDeltaPPM(otherInTwoPDic)
    # writeDeltaPPM(xlink2deltaDic, 'xlink2ppm.txt')
    # writeDeltaPPM(mono2deltaDic, 'mono2ppm.txt')
    # writeDeltaPPM(other2deltaDic, "other2ppm.txt")
    # xlinkInOnePdic, monoInOnePdic, otherInOnePdic = stacOnePairinfo(ms2TriggerMS3dic)
    # print(len(xlinkInOnePdic), len(monoInOnePdic), len(otherInOnePdic))
    # xlink1deltaDic = calPairedIonDeltaPPM(xlinkInOnePdic)
    # mono1deltaDic = calPairedIonDeltaPPM(monoInOnePdic)
    # other1deltaDic = calPairedIonDeltaPPM(otherInOnePdic)
    # writeDeltaPPM(xlink1deltaDic, "xlink1ppm.txt")
    # writeDeltaPPM(mono1deltaDic, "mono1ppm.txt")
    # writeDeltaPPM(other1deltaDic, "other1ppm.txt")


if __name__ == "__main__":
    main()