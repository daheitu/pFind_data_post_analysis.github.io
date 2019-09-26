# coding = utf-8

import os, sys, time
import copy
from numpy import zeros
sys.path.append(
    r"F:\OneDrive\github\pFind_data_post_analysis.github.io\dsso_feature_finder"
)
from filter_Reuslt import main

os.chdir(r"F:\MS_DATA_STORAGE\20190819\multiNCE\inclu_random_reverse")

flPath = r"./BSA_DSSO_INCLU_RANDOM_REVERSE_R1.ms2"
incluPath = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.08.19.19.57.19\inclusion_list.csv"
mgfPath = r"./BSA_DSSO_INCLU_RANDOM_REVERSE_R1.mgf"

LinkerMass = 158.004  # 交联剂质量
longArmMass = 85.9826
shortArmMass = 54.0106

mpMassTable={'A':71.037114, 'R':156.101111, 'N':114.042927, 'D':115.026943,\
    'C':103.009185, 'E':129.042593,'Q':128.058578,'G':57.021464, 'H':137.058912,\
    'I':113.084064, 'L':113.084064,'K':128.094963,'M':131.040485, 'F':147.068414,\
    'P':97.052764, 'S':87.032028, 'T':101.047679, 'U':150.95363, 'W':186.079313,\
    'Y':163.06332, 'V':99.068414, 'H2O':18.01056,'H1':1.00782}

mpModMass = {
    'Carbamidomethyl[C]': 57.021464,
    'Oxidation[M]': 15.994915
}  # 添加的修饰名称及质量
mstol = 20.0


def calMassPepPlusMod(massList, modMass):
    pepMass = sum(massList) + mpMassTable['H2O']
    return pepMass + modMass + mpMassTable['H1']


# 根据肽段和修饰生成理论碎片离子
def getTheroIons(pep, modCell):
    a, b = pep.split('-')

    site_a = int(''.join([i for i in a if i.isdigit()])) - 1
    site_b = int(''.join([i for i in b if i.isdigit()])) - 1

    pep_a = a.split("(")[0]  #''.join([i for i in a if i.isalpha()])
    pep_b = b.split("(")[0]  # ''.join([i for i in b if i.isalpha()])

    pep_a_len = len(pep_a)
    pep_b_len = len(pep_b)

    a_mod_mass_list = [0] * pep_a_len
    b_mod_mass_list = [0] * pep_b_len

    if modCell != 'null':
        modList = modCell.split(';')
        for m in modList:
            name, pos = m[:-1].split('(')
            pos = int(pos)
            if pos <= pep_a_len + 1:
                pos -= 1  # 与python的index 从0开始对齐
                if pos == pep_a_len:  #考虑修饰在C端的情况此时pos = peptide——A长度
                    pos -= 1  #超出了list的index，因此减一
                elif pos == -1:  # 考虑修饰在N端的情况 pos = -1，因此+1，算到第一个氨基酸是
                    pos += 1
                a_mod_mass_list[pos] += mpModMass[name]
            else:
                pos = pos - pep_a_len - 4
                if pos == -1:
                    pos = 0
                elif pos == pep_b_len:
                    pos -= 1
                b_mod_mass_list[pos] += mpModMass[name]

    a_mass_list = [0] * pep_a_len  # 残基质量
    b_mass_list = [0] * pep_b_len

    for i in range(pep_a_len):
        a_mass_list[i] = mpMassTable[a[i]] + a_mod_mass_list[i]
    for i in range(pep_b_len):
        b_mass_list[i] = mpMassTable[b[i]] + b_mod_mass_list[i]

    a_long_mass = calMassPepPlusMod(a_mass_list, longArmMass)

    a_short_mass = calMassPepPlusMod(a_mass_list, shortArmMass)

    b_long_mass = calMassPepPlusMod(b_mass_list, longArmMass)

    b_short_mass = calMassPepPlusMod(b_mass_list, shortArmMass)

    return [a_long_mass, a_short_mass, b_long_mass, b_short_mass]


# 读取xlpep文件，提取所有ID谱图
def getIDspec(xlpepPath):
    idSpecList = []
    f = open(xlpepPath, 'r').readlines()
    for line in f[2:]:
        lineList = line.split(",")
        if lineList[0] == "" and lineList[1].isdigit():
            idSpecList.append(lineList[2])
    idSpecList = sorted(idSpecList, key = lambda x: int(x.split(".")[1]))
    return idSpecList


def getIDspecMZdic(idSpecList, mgfPath):
    IDspecMZdic = {}
    endidListScan = int(idSpecList[-1]split(".")[1])
    mgf = open(mgfPath, 'r').readlines()
    i = 0
    while i < len(mgf):
        if mgf[i].starswith("TITLE"):
            title = mgf[i][6:-1]
            scan = int(title.split(".")[1])
            if scan > endidListScan:
                break
            else:
                if title not in idSpecList:
                    i += 1
                else:
                    p = i + 4
                    ms2Dic = {}
                    while p < len(mgf):
                        if mgf[p].starswith("END IONS"):
                            break
                        else:
                            mz, ints = mgf[:-1].split(" ")
                            ms2Dic[mz] = ints
                    IDspecMZdic[title] = ms2Dic
                    i = p + 1
        else:
            i += 1            
    return IDspecMZdic


def isMatched(mz, spec):
    lowTgtMZ, upTgtMZ = generate_ion_mass_range(mz)
    for i in range(len(spec)):
        if 'END' in spec[i]:
            break
        exp_mz = float(spec[i].split(' ')[0])
        if exp_mz > upTgtMZ:
            break
        elif exp_mz >= lowTgtMZ:
            return True
        else:
            continue

    return False


#生成给定mz的容忍的上下质量区间
def generate_ion_mass_range(num):
    deta = num * mstol / 1000000
    return num - deta, num + deta


# 对给定MZ。计算是否在谱图中出现，如果出现，返回匹配值和其强度
def isMatchForReport(mz, ms2_dic):
    mz_list = sorted(list(ms2_dic.keys()))
    max_ins = max(list(ms2_dic.values()))
    lowTgtMZ, upTgtMZ = generate_ion_mass_range(mz)
    if lowTgtMZ > mz_list[-1]:
        return False
    else:
        i = 0
        while i < len(mz_list):
            if mz_list[i] < lowTgtMZ:
                i += 1
            elif mz_list[i] <= upTgtMZ:
                return True, mz_list[i], ms2_dic[mz_list[i]] / max_ins
            else:
                return False


#对给定的报告离子的1+质量列表，计算多价态下其是否存在，并写入匹配文件
def match_report_ion(mass_list, spec, max_c, fmatched):
    ms2_dic = {}
    for i in range(len(spec) - 1):
        ms2_mz = float(spec[i][:-1].split(" ")[0])
        ms2_ins = float(spec[i][:-1].split(" ")[1])
        ms2_dic[ms2_mz] = ms2_ins
    rep_mask = [0] * 4
    match_dic = {}
    mz_list = sorted(list(ms2_dic.keys()))
    max_ins = max(list(ms2_dic.values()))
    for c in range(1, max_c + 1):
        mass_cur_list = [(m + 1.0078 * (c - 1)) / c for m in mass_list]
        for i in range(len(mass_cur_list)):
            if isMatchForReport(mass_cur_list[i], ms2_dic):
                matchedMZ, matchedRelInts = isMatchForReport(
                    mass_cur_list[i], ms2_dic)[1:]
                rep_mask[i] += 1
                if i == 0:
                    fmatched.write('charge=%d\tmz=%f\ttype=Alpha_long_mass\n' %
                                   (c, mass_list[i]))
                elif i == 1:
                    fmatched.write(
                        'charge=%d\tmz=%f\ttype=Alpha_short_mass\n' %
                        (c, mass_list[i]))
                elif i == 2:
                    fmatched.write('charge=%d\tmz=%f\ttype=Beta_long_mass\n' %
                                   (c, mass_list[i]))
                else:
                    fmatched.write('charge=%d\tmz=%f\ttype=Beta_short_mass\n' %
                                   (c, mass_list[i]))

                if i not in match_dic:
                    match_dic[i] = [[
                        c, mass_cur_list[i], matchedMZ, matchedRelInts
                    ]]
                else:
                    match_dic[i].append(
                        [c, mass_cur_list[i], matchedMZ, matchedRelInts])
    return match_dic


#对于报告离子的matchdic，检测是否成对
def detectPAIRreporter(match_dic, m, n):
    pairRepIon_list = []
    if m in match_dic and n in match_dic:
        chrgm_list = [ele[0] for ele in match_dic[m]]
        chrgn_list = [ele[0] for ele in match_dic[n]]
        intsm_list = [ele[-1] for ele in match_dic[m]]
        intsn_list = [ele[-1] for ele in match_dic[n]]

        if len(chrgm_list) <= len(chrgn_list):
            for chg in chrgm_list:
                if chg in chrgn_list:
                    ave_ints = (intsm_list[chrgm_list.index(chg)] +
                                intsn_list[chrgn_list.index(chg)]) / 2
                    pairRepIon_list.append([chg, ave_ints])
        else:
            for chg in chrgn_list:
                if chg in chrgm_list:
                    ave_ints = (intsm_list[chrgm_list.index(chg)] +
                                intsn_list[chrgn_list.index(chg)]) / 2
                    pairRepIon_list.append([chg, ave_ints])

    if pairRepIon_list:
        pairRepIon_list = sorted(pairRepIon_list, key=lambda x: x[-1])

    return pairRepIon_list


def summary_rep_macthDic(match_dic):
    pairAlphaRep_list = detectPAIRreporter(match_dic, 0, 1)
    pairBetaRep_list = detectPAIRreporter(match_dic, 2, 3)

    xAlphapair_bool = False
    xBetaPair_bool = False
    maxBetaRepInts = "Na"
    maxAlphaRepInts = "Na"
    if pairAlphaRep_list:
        xAlphapair_bool = True
        maxAlphaRepInts = pairAlphaRep_list[-1][-1]

    if pairBetaRep_list:
        xBetaPair_bool = True
        maxBetaRepInts = pairBetaRep_list[-1][-1]

    pair_num = [xBetaPair_bool, xAlphapair_bool].count(True)
    if pair_num > 0:
        detect_bool = True
    else:
        detect_bool = False

    return detect_bool, pair_num, maxBetaRepInts, maxAlphaRepInts


def calIonRatio(pep, spec, charge, modCell, fmatched):
    max_c = charge -1 
    a_long_mass, a_short_mass, b_long_mass, b_short_mass = getTheroIons(pep, modCell)
    reporter_ions= [a_long_mass, a_short_mass, b_long_mass, b_short_mass]

    fmatched.write('---------reporter ions---------------\n')
    match_dic = match_report_ion(reporter_ions, spec, max_c, fmatched)
    detect_bool, pair_num, maxBetaRepInts, maxAlphaRepInts = summary_rep_macthDic(
        match_dic)

    return  detect_bool, pair_num, maxBetaRepInts, maxAlphaRepInts, 


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
                            return mzdic[mzList[i], charge]
                        else:
                            return False
                    else:
                        i += 1


def getScanNCE(ms2flpath):
    scanNCEdic = {}
    f = open(ms2flpath, 'r').readlines()
    scanNCEdic = {}
    i = 0
    while i < len(f):
        if f[i][0] != "S":
            i += 1
        else:
            scan = int(f[i].split("\t")[1])
            mz = float(f[i].split("\t")[-1])
            charge = f[i + 9].split("\t")[1]
            nce = f[i + 6].split("\t")[2].split(" ")[7].split("@")[1][3:5]
            scanNCEdic[scan] = [mz, charge, nce]
            i += 10
    return scanNCEdic


def getmzPepMod(incluPath):
    mzPepModDic = {}
    f = open(incluPath, 'r').readlines()
    for line in f[1:]:
        lineList = line.split(',')
        charge = lineList[1]
        therMZ = float(lineList[3])
        pep = lineList[0]
        mod = lineList[2]
        #print(therMZ, charge)
        mzPepModDic[therMZ, charge] = [pep, mod]
    return mzPepModDic


def pepModToTheoIon(mzPepModDic, ftheo_ions):
    pepModToTheoIonDic = {}
    pepModList = []
    for key in mzPepModDic:
        [pep, mod] = mzPepModDic[key]
        if (pep, mod) not in pepModToTheoIonDic:
            by_ions, reporter_ions, intXPions = getTheroIons(pep, mod)
            pepModToTheoIonDic[pep, mod] = [by_ions, reporter_ions, intXPions]
            writeTheorIon(pep, mod, by_ions, reporter_ions, intXPions,
                          ftheo_ions)

    return pepModToTheoIonDic


def linkscanNCEmzPepMod():
    scanNCEdic = getScanNCE(flPath)
    mzPepModDic = getmzPepMod(incluPath)
    delScans = []
    for scan in scanNCEdic:
        realMZ, charge = scanNCEdic[scan][:2]
        judgeBool = findMZinmzDic(realMZ, charge, mzPepModDic)
        if judgeBool:
            scanNCEdic[scan].extend(judgeBool)
        else:
            print(scan, realMZ, charge)
            delScans.append(scan)
    for scan in delScans:
        del scanNCEdic[scan]
    return scanNCEdic


def main():
    bg_time = time.localtime(time.time())

    fmatched = open("./matched_info", 'w')
    ftheo_ions = open("./theoretical_ions", 'w')
    mzPepModDic = getmzPepMod(incluPath)
    pepModToTheoIonDic = pepModToTheoIon(mzPepModDic, ftheo_ions)
    scanNCEdic = linkscanNCEmzPepMod()

    repDic = {}
    scanSpecdic = getTargetSpec(mgfPath)
    for scan in scanNCEdic:
        ttl = scanNCEdic[scan]
        print(scan)
        pep = ttl[3]
        modCell = ttl[4]
        charge = int(ttl[1])
        nce = ttl[2]

        spec = scanSpecdic[scan]

        extendInfo = calIonRatio(pep, spec, charge, modCell,
                                 pepModToTheoIonDic, fmatched)
        ttl.extend(extendInfo)
        ttl.insert(0, scan)
        repDic[scan] = ttl

    fmatched.close()
    ftheo_ions.close()

    b = open("reportFile_commBY.csv", 'w')
    b.write("scan, mz, charge, nce, pep, mods, regular_by_ion_ratio,\
        xlink_by_ion_ratio, detect_bool, pair_num, rep_ints_beta,\
            rep_ints_alpha, intPXratio, commBYalpha, commBYbeta" + "\n")

    for scan in repDic:
        wlist = [str(ele) for ele in repDic[scan]]
        b.write(",".join(wlist) + "\n")
    b.close()

    endTime = time.localtime(time.time())
    print(bg_time, "\n", endTime)


if __name__ == '__main__':
    main()
#main()