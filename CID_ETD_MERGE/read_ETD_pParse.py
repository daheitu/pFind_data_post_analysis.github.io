import os

os.chdir(r"G:\DSSO_1008\CID_BASED_Methods\ETD")


def averageTwoNum(num1, num2):
    # digitLen = len(str(num1).split('.')[1])
    digitLen = 5
    return round((num1+num2)/2, digitLen)


def generate_ion_mass_range(num, tol_ppm):
    deta = num * tol_ppm / 1000000
    return num - deta, num + deta


def compareTwoMZ(mz1, mz2, tol_ppm):
    lowMZ1,upMZ1 = generate_ion_mass_range(mz1, tol_ppm)
    if mz2 > lowMZ1 and mz2 < upMZ1:
        return True
    else:
        return False


def detailCompTwoMZ(mz1, mz2, tol_ppm):
    if mz1 > mz2:
        if compareTwoMZ(mz1, mz2, tol_ppm):
            return 0
        else:
            return -1
    else:
        if compareTwoMZ(mz1, mz2, tol_ppm):
            return 0
        else:
            return 1


def readETDinfo(path):
    etd = open(path, 'r').readlines()
    info_dic = {}
    i = 0
    while i < len(etd):
        if etd[i].strip() == "BEGIN IONS":
            title = etd[i+1].strip().split("=")[1]
            scanNum = int(title.split('.')[1])
            rawName = title.split('.')[0]
            charge = int("".join([i for i in etd[i+2] if i.isdigit()]))
            rtLine = etd[i+3]
            preMZ = etd[i+4].strip().split("=")[1]
            info_dic[scanNum] = [charge, preMZ]
        else:
            print('wrong')
        ms2MzIntDic = {}
        p = i+5
        while p < len(etd):
            if etd[p].strip() == "END IONS":
                break
            else:
                if etd[p][0].isdigit():
                    [mz, intes] = etd[p].strip().split(" ")
                    ms2MzIntDic[float(mz)] = float(intes)
                else:
                    pass
            p += 1
        
        info_dic[scanNum] = [charge, preMZ, ms2MzIntDic, rtLine, rawName]
        i = p + 1
    return info_dic

# print(info_dic)

def readCIDinfo(path):
    cid = open(path, 'r').readlines()
    info_dic = {}
    i = 0
    while i < len(cid):
        if cid[i].strip() == "BEGIN IONS":
            title = cid[i+1].strip().split("=")[1]
            scanNum = int(title.split('.')[1])
            charge = int("".join([i for i in cid[i+2] if i.isdigit()]))
            mz = cid[i+4].strip().split("=")[1]
            info_dic[scanNum] = [charge, mz]
        else:
            print('wrong')
        p = i+5
        while p < len(cid):
            if cid[p].strip() == "END IONS":
                break
            else:
                pass 
            p += 1
        i = p + 1
    return info_dic


def writePreInfo(chrg, mz, rtLine, scanETD, scanCID, rawName, openedFile):
    openedFile.write('BEGIN IONS' + '\n')
    titleList = ['TITLE='+rawName, str(scanCID), str(scanETD), str(chrg), "0.dta"]
    openedFile.write('.'.join(titleList) + '\n')
    chrgList = ['CHARGE=', str(chrg), '+']
    openedFile.write(''.join(chrgList)+'\n')
    openedFile.write(rtLine)
    openedFile.write("PEPMASS=" + mz + '\n')


def writeMergedInfo(etdMS2dic, cidMS2dic, openedFile):
    mergedMs2Dic = {}
    etdMZlist = list(etdMS2dic.keys())
    maxETDIntes = max(list(etdMS2dic.values()))
    cidMZlist = list(cidMS2dic.keys())
    maxCIDIntes = max(list(etdMS2dic.values()))
    averMaxIntes = averageTwoNum(maxETDIntes, maxCIDIntes)
    etdMZlist.sort()
    cidMZlist.sort()
    m = 0; n = 0
    while m < len(etdMZlist) and n < len(cidMZlist):
        mzETD = etdMZlist[m]; mzCID = cidMZlist[n]
        compRe = detailCompTwoMZ(mzETD, mzCID, 10) 
        if compRe == 0:
            m += 1; n += 1
            averMZ = averageTwoNum(mzETD, mzCID)
            addIntes = averageTwoNum(etdMS2dic[mzETD], cidMS2dic[mzCID])*2
            mergedMs2Dic[averMZ] = addIntes
            openedFile.write(str(averMZ) + " " +str(addIntes) + "\n")
        elif compRe == -1:
            n += 1
            mergedMs2Dic[mzCID] = cidMS2dic[mzCID]
            caliIntes = round(cidMS2dic[mzCID] / maxCIDIntes * averMaxIntes, 1) 
            openedFile.write(' '.join([str(mzCID),  str(caliIntes)]) + '\n')
        else:
            m += 1
            mergedMs2Dic[mzETD] = etdMS2dic[mzETD]
            caliIntes = round(etdMS2dic[mzETD] / maxETDIntes * averMaxIntes, 1)
            openedFile.write(str(mzETD) + ' ' + str(caliIntes) + '\n')
    mergedMZlist = list(mergedMs2Dic.keys())
    #print(etdMZlist)
    #print(cidMZlist)
    #print(mergedMZlist)
    openedFile.write("END IONS" + '\n')
    

def main():
    repList=[]
    delta_num = []
    etdDic = readETDinfo("DSSO_CID_ETD_MS2_RE23_MIPS_T1_ETDFT.mgf")
    cidDic = readETDinfo("DSSO_CID_ETD_MS2_RE23_MIPS_T1_CIDFT.mgf")
    mrg = open("DSSO_CID_ETD_MS2_RE23_MIPS_T1_CIDETDFT.mgf", 'w')
    
    for num in etdDic:
        matchETD2CIDbool = False
        etdCharge, etdMZ, etdMzIntesDic, rtLine, rawName = etdDic[num]
        for k in range(1, 15):
            if num - k not in cidDic:
                pass
            else:
                cidCharge, cidMZ, cidMzIntesDic = cidDic[num-k][:3]
                compPrec = compareTwoMZ(float(etdMZ), float(cidMZ), 20)
                if etdCharge == cidCharge and compPrec:
                    matchETD2CIDbool = True
                    pair = [num, num - k]
                    repList.append(pair)
                    delta_num.append(k)
                    break
                else:
                    pass
        if matchETD2CIDbool == True:
            writePreInfo(etdCharge, etdMZ, rtLine, num, num-k, rawName, mrg)
            writeMergedInfo(etdMzIntesDic, cidMzIntesDic, mrg)
    
    mrg.close()
    delta_num.sort()
    # print(repList)
    print(max(delta_num))
    print(len(etdDic))
    print(len(delta_num)/len(etdDic))


if __name__ == "__main__":
    main()
