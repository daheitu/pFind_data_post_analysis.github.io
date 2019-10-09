# coding = utf-8

import os, sys, time
import copy
from numpy import zeros
sys.path.append(
    r"F:\OneDrive\github\pFind_data_post_analysis.github.io\dsso_feature_finder"
)
from filter_Reuslt import main

#os.chdir(r"F:\MS_DATA_STORAGE\20190819\multiNCE\inclu_random_reverse")

flPath = r"./BSA_DSSO_INCLU_RANDOM_REVERSE_R1.ms2"
incluPath = r"./inclusion_list.csv"
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


def calPepBYions(massList, linkSite, linkAddMass):
    pepMass = sum(massList) + mpMassTable['H2O'] + linkAddMass\
         + 2 * mpMassTable['H1']
    pep_len = len(massList)
    tmp_mass = mpMassTable['H1']
    b1Ion_pep = [0] * (pep_len - 1)
    y1Ion_pep = [0] * (pep_len - 1)
    for i in range(pep_len - 1):
        if i == linkSite:
            tmp_mass += linkAddMass
        else:
            pass
        tmp_mass += massList[i]
        b1Ion_pep[i] = tmp_mass
        y1Ion_pep[i] = pepMass - tmp_mass

    return b1Ion_pep, y1Ion_pep


def calInternalIon(massList, linkSite, linkAddMass):
    pep_len = len(massList)
    Intenal_Pep_List = [0] * (pep_len - 1)
    for i in range(pep_len - 1):
        if i < linkSite:
            Intenal_Pep_List[i] = sum(massList[i + 1:]) + mpMassTable['H2O']\
                 + linkAddMass + 1.00782
        else:
            Intenal_Pep_List[i] = sum(massList[:i + 1]) + linkAddMass + 1.00782
    return Intenal_Pep_List


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

    ################### 4 reporter ions mass 1+ ############################
    a_long_mass = calMassPepPlusMod(a_mass_list, longArmMass)

    a_short_mass = calMassPepPlusMod(a_mass_list, shortArmMass)

    b_long_mass = calMassPepPlusMod(b_mass_list, longArmMass)

    b_short_mass = calMassPepPlusMod(b_mass_list, shortArmMass)
    ###################4 reporter ions mass 1+ ############################

    ######################regular b,y ions 1+########################
    pep_mass = a_long_mass + b_short_mass + mpMassTable['H2O']
    a_linker_mass = sum(a_mass_list) + mpMassTable['H2O'] + LinkerMass
    b_linker_mass = sum(b_mass_list) + mpMassTable['H2O'] + LinkerMass

    b1IonPepAlpha, y1IonPepAlpha = calPepBYions(a_mass_list, site_a,
                                                b_linker_mass)

    b1IonPepBeta, y1IonPepBeta = calPepBYions(b_mass_list, site_b,
                                              a_linker_mass)
    ######################regular b,y ions 1+########################

    #############Internal ions mass 1+  cleave at linker and backbone######
    IntenalIonLongAlpha = calInternalIon(a_mass_list, site_a, longArmMass)
    IntenalIonShortAlpha = calInternalIon(a_mass_list, site_a, shortArmMass)

    IntenalIonLongBeta = calInternalIon(b_mass_list, site_b, longArmMass)
    IntenalIonShortBeta = calInternalIon(b_mass_list, site_b, shortArmMass)
    #############Internal ions mass 1+  cleave at linker and backbone######

    return [b1IonPepAlpha, y1IonPepAlpha, b1IonPepBeta, y1IonPepBeta], \
            [a_short_mass, a_long_mass, b_short_mass, b_long_mass],\
            [IntenalIonLongAlpha, IntenalIonShortAlpha, IntenalIonLongBeta, IntenalIonShortBeta]


def getTargetSpec(spec_path):
    f = open(spec_path).readlines()
    scanSpecdic = {}
    i = 0
    while i < len(f):
        if f[i].strip() != "BEGIN IONS":
            print(i)
        else:
            scan = int(f[i + 2].split("=")[1].strip())
            vOneSpec = []
            p = i + 6
            while p < len(f):
                if f[p].strip() == "END IONS":
                    vOneSpec.append(f[p])
                    break
                else:
                    vOneSpec.append(f[p])
                    p += 1
            scanSpecdic[scan] = vOneSpec
            i = p + 2
    return scanSpecdic


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
                    fmatched.write('charge=%d\tmz=%f\ttype=Alpha_short_mass\t1+mass=%f\n' %
                                   (c, matchedMZ, mass_list[i]))
                elif i == 1:
                    fmatched.write(
                        'charge=%d\tmz=%f\ttype=Alpha_long_mass\t1+mass=%f\n' %
                        (c, matchedMZ, mass_list[i]))
                elif i == 2:
                    fmatched.write('charge=%d\tmz=%f\ttype=Beta_short_mass\t1+mass=%f\n' %
                                   (c, matchedMZ, mass_list[i]))
                else:
                    fmatched.write('charge=%d\tmz=%f\ttype=Beta_long_mass\t1+mass=%f\n' %
                                   (c, matchedMZ, mass_list[i]))

                if i not in match_dic:
                    match_dic[i] = [[
                        c, mass_cur_list[i], matchedMZ, matchedRelInts
                    ]]
                else:
                    match_dic[i].append(
                        [c, mass_cur_list[i], matchedMZ, matchedRelInts])
    fmatched.write("reporter ions matched: " + str(rep_mask) + "\n")
    return match_dic, rep_mask


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


#汇总报告离子
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


def matchBYions(b1, y1, site, max_c, spec, fmatched):
    common_mask = [0] * len(b1)
    xlink_mask = [0] * len(y1)

    for c in range(1, max_c + 1):
        b_cur = []
        y_cur = []

        b_cur = [(m + 1.0078 * (c - 1)) / float(c) for m in b1]
        y_cur = [(m + 1.0078 * (c - 1)) / float(c) for m in y1]

        for i in range(len(b1)):
            if i < site and isMatched(b_cur[i], spec):
                common_mask[i] = 1
                fmatched.write('site=%d\tcharge=%d\tmz=%f\ttype=common_b\n' %
                               (i, c, b_cur[i]))
            if i >= site and isMatched(y_cur[i], spec):
                common_mask[i] = 1
                fmatched.write('site=%d\tcharge=%d\tmz=%f\ttype=common_y\n' %
                               (i, c, y_cur[i]))
            if i < site and isMatched(y_cur[i], spec):
                xlink_mask[i] = 1
                fmatched.write('site=%d\tcharge=%d\tmz=%f\ttype=xlink_y\n' %
                               (i, c, y_cur[i]))
            if i >= site and isMatched(b_cur[i], spec):
                xlink_mask[i] = 1
                fmatched.write('site=%d\tcharge=%d\tmz=%f\ttype=xlink_b\n' %
                               (i, c, b_cur[i]))
    fmatched.write('common ion matched: ' + str(common_mask) + '\n')
    fmatched.write('xlink ion matched: ' + str(xlink_mask) + '\n')
    return common_mask, xlink_mask


def judgeInternalPair(inLongMusk, inShortMusk):
    internalPAIRmusk = [0] * len(inLongMusk)
    if len(inLongMusk) != len(inShortMusk):
        print("Warning: inLongMusk not equel to inShortMusk")
    else:
        for i in range(len(inLongMusk)):
            if inLongMusk[i] > 0 and inShortMusk[i] > 0:
                internalPAIRmusk[i] = 1
            else:
                internalPAIRmusk[i] = 0
    return internalPAIRmusk


def matchPepIntXP(In_pepA_long, In_pepA_short, In_pepB_long, In_pepB_short,
                  spec, fmatched):
    max_c = 2
    inPAXL = copy.deepcopy(In_pepA_long)
    inPAXS = copy.deepcopy(In_pepA_short)
    inPBXL = copy.deepcopy(In_pepB_long)
    inPBXS = copy.deepcopy(In_pepB_short)
    inPAXL_mask = [0] * len(inPAXL)
    inPAXS_mask = [0] * len(inPAXS)
    inPBXL_mask = [0] * len(inPBXL)
    inPBXS_mask = [0] * len(inPBXS)
    pairAlphaMusk = [0] * len(inPAXL)
    pairBetaMusk = [0] * len(inPBXL)

    for c in range(1, max_c + 1):
        inPAXL_cur = [(m + 1.0078 * (c - 1)) / c for m in inPAXL]
        inPAXS_cur = [(m + 1.0078 * (c - 1)) / c for m in inPAXS]
        inPBXL_cur = [(m + 1.0078 * (c - 1)) / c for m in inPBXL]
        inPBXS_cur = [(m + 1.0078 * (c - 1)) / c for m in inPBXS]

        for i in range(len(inPAXL_cur)):
            if isMatched(inPAXL_cur[i], spec):
                inPAXL_mask[i] += 1
                fmatched.write(
                    'site=(%d, AL)\tcharge=%d\tmz=%f\ttype=In_pepA_long\n' %
                    (i, c, inPAXL_cur[i]))

            if isMatched(inPAXS_cur[i], spec):
                inPAXS_mask[i] += 1
                fmatched.write(
                    'site=(%d, AS)\tcharge=%d\tmz=%f\ttype=In_pepA_short\n' %
                    (i, c, inPAXS_cur[i]))
            
            if isMatched(inPAXL_cur[i], spec) and isMatched(inPAXS_cur[i], spec):
                pairAlphaMusk[i] += 1

        for i in range(len(inPBXL_cur)):
            if isMatched(inPBXL_cur[i], spec):
                inPBXL_mask[i] += 1
                fmatched.write(
                    'site=(%d, BL)\tcharge=%d\tmz=%f\ttype=In_pepB_long\n' %
                    (i, c, inPBXL_cur[i]))
            if isMatched(inPBXS_cur[i], spec):
                inPBXS_mask[i] += 1
                fmatched.write(
                    'site=(%d, BS)\tcharge=%d\tmz=%f\ttype=In_pepB_short\n' %
                    (i, c, inPBXS_cur[i]))
            
            if isMatched(inPBXL_cur[i], spec) and isMatched(inPBXS_cur[i], spec):
                pairBetaMusk[i] += 1
    
    fmatched.write("Alpha interPX Long match: " + str(inPAXL_mask) + "\n")
    fmatched.write("Alpha interPX short match: " + str(inPAXS_mask) + "\n")
    fmatched.write("Beta interPX Long match: " + str(inPBXL_mask) + "\n")
    fmatched.write("Beta interPX short match: "+ str(inPBXS_mask) + "\n")
    alphaLen = len(In_pepA_long)
    betaLen = len(In_pepB_long)
    longCoverRatio = 1 - (inPAXL_mask.count(0) + inPBXL_mask.count(0))/ (alphaLen + betaLen)
    shortCoverRatio = 1 - (inPAXS_mask.count(0) + inPBXS_mask.count(0)) / (alphaLen + betaLen)
    total_unmat = inPAXL_mask.count(0) + inPAXS_mask.count(
        0) + inPBXL_mask.count(0) + inPBXS_mask.count(0)
    intXP_ratio = 1 - total_unmat / float(len(inPAXL) + len(inPBXL)) / 2
    pairRATIO = 1- (pairAlphaMusk.count(0) + pairBetaMusk.count(0))/ (alphaLen + betaLen)
    return intXP_ratio, longCoverRatio, shortCoverRatio, pairRATIO  


def writeTheorIon(pep, mod, by_ions, reporter_ions, intXPions, ftheo_ions):
    b1_a, y1_a, b1_b, y1_b = by_ions

    a_long_mass, a_short_mass, b_long_mass, b_short_mass = reporter_ions

    In_pepA_long, In_pepA_short, In_pepB_long, In_pepB_short = intXPions

    ftheo_ions.write('XLPeptide=%s\tmodification=%s\n' % (pep, mod))
    ftheo_ions.write('--------alpha backbone ions (alpha b1+)---------\n')
    ftheo_ions.write("\t".join([str(ele)
                                for ele in b1_a if type(ele) != str]) + "\n")

    ftheo_ions.write('--------alpha backbone ions (alpha y1+)---------\n')
    ftheo_ions.write("\t".join([str(ele)
                                for ele in y1_a if type(ele) != str]) + "\n")

    ftheo_ions.write('--------beta backbone ions (beta b1+)---------\n')
    ftheo_ions.write("\t".join([str(ele)
                                for ele in b1_b if type(ele) != str]) + "\n")

    ftheo_ions.write('--------beta backbone ions (beta y1+)---------\n')
    ftheo_ions.write("\t".join([str(ele)
                                for ele in y1_b if type(ele) != str]) + "\n")

    ftheo_ions.write('--------reporter ions (alpha long 1+)---------\n')
    ftheo_ions.write(str(a_long_mass) + "\n")

    ftheo_ions.write('--------reporter ions (alpha short 1+)---------\n')
    ftheo_ions.write(str(a_short_mass) + "\n")

    ftheo_ions.write('--------reporter ions (beta long 1+)---------\n')
    ftheo_ions.write(str(b_long_mass) + "\n")

    ftheo_ions.write('--------reporter ions (beta short 1+)---------\n')
    ftheo_ions.write(str(b_short_mass) + "\n")


def calIonRatio(pep, spec, charge, modCell, pepModToTheoIonDic, fmatched):
    a, b = pep.split('-')

    site_a = int(''.join([i for i in a if i.isdigit()])) - 1
    site_b = int(''.join([i for i in b if i.isdigit()])) - 1

    a = ''.join([i for i in a if i.isalpha()])
    b = ''.join([i for i in b if i.isalpha()])

    pep_a_len = len(a)
    pep_b_len = len(b)

    max_c = min(3, charge)

    by_ions, reporter_ions, intXPions = pepModToTheoIonDic[(pep, modCell)]

    b1_a, y1_a, b1_b, y1_b = by_ions
    a_long_mass, a_short_mass, b_long_mass, b_short_mass = reporter_ions

    In_pepA_long, In_pepA_short, In_pepB_long, In_pepB_short = intXPions

    fmatched.write('XLpep=%s\tmod=%s\tcharge=%d\n' % (pep, modCell, charge))
    fmatched.write('---------alpha---------------\n')
    a_common_mask, a_xlink_mask = matchBYions(b1_a, y1_a, site_a, max_c, spec,
                                              fmatched)

    fmatched.write('---------beta---------------\n')
    b_common_mask, b_xlink_mask = matchBYions(b1_b, y1_b, site_b, max_c, spec,
                                              fmatched)

    fmatched.write('---------reporter ions---------------\n')
    match_dic, rep_mask = match_report_ion(reporter_ions, spec, max_c, fmatched)
    detect_bool, pair_num, maxBetaRepInts, maxAlphaRepInts = summary_rep_macthDic(
        match_dic)
    repState = ";".join([str(ele) for ele in rep_mask])

    fmatched.write('---------internal PX Ions---------------\n')
    intXP_ratio, longCoverRatio, shortCoverRatio, pairRATIO = matchPepIntXP(In_pepA_long, In_pepA_short, In_pepB_long,
                               In_pepB_short, spec, fmatched)

    sumCommAlpha = sum(a_common_mask)
    sumXlinkAlpha = sum(a_xlink_mask)
    sumCommBeta = sum(b_common_mask)
    sumXlinkBeta = sum(b_xlink_mask)

    regular_by_ion_ratio = (sumCommAlpha + sumCommBeta) / (pep_b_len +
                                                           pep_a_len - 2)
    xlink_by_ion_ratio = (sumXlinkAlpha + sumXlinkBeta) / (pep_b_len +
                                                           pep_a_len - 2)
    commByAlpha = sumCommAlpha / (pep_a_len - 1)
    commByBeta = sumCommBeta / (pep_b_len - 1)

    return regular_by_ion_ratio, xlink_by_ion_ratio, commByAlpha, commByBeta,\
        detect_bool, pair_num, repState, maxBetaRepInts, maxAlphaRepInts,\
        intXP_ratio, longCoverRatio, shortCoverRatio, pairRATIO


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
        fmatched.write("scan number= %d" % scan + "\n")
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

    b = open("reportFile_commBY_v2.csv", 'w')
    b.write("scan, mz, charge, nce, pep, mods, regular_by_ion_ratio, \
        xlink_by_ion_ratio, commByAlpha, commByBeta, detect_bool, \
        pair_num, repState, maxBetaRepInts, maxAlphaRepInts,\
        intXP_ratio, longCoverRatio, shortCoverRatio, pairRATIO" + "\n")

    for scan in repDic:
        wlist = [str(ele) for ele in repDic[scan]]
        b.write(",".join(wlist) + "\n")
    b.close()

    endTime = time.localtime(time.time())
    print(bg_time, "\n", endTime)


if __name__ == '__main__':
    main()
#main()