# coding = utf-8

import os, sys
import copy
from numpy import zeros
sys.path.append(r"F:\OneDrive\github\pFind_data_post_analysis.github.io\dsso_feature_finder")
from filter_Reuslt import main

os.chdir(r"F:\MS_DATA_STORAGE\20190819\multiNCE\inclu_random_reverse")

flPath = r"./BSA_DSSO_INCLU_RANDOM_REVERSE_R1.ms2"
incluPath = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.08.19.19.57.19\inclusion_list.csv"
mgfPath = r"./BSA_DSSO_INCLU_RANDOM_REVERSE_R1.mgf"
fmatched = open("./matched_info", 'w')
ftheo_ions = open("./theoretical_ions", 'w')

LinkerMass = 158.004  # 交联剂质量
longArmMass = 85.9826
shortArmMass = 54.0106

mpMassTable={'A':71.037114, 'R':156.101111, 'N':114.042927, 'D':115.026943,\
    'C':103.009185, 'E':129.042593,'Q':128.058578,'G':57.021464, 'H':137.058912,\
    'I':113.084064, 'L':113.084064,'K':128.094963,'M':131.040485, 'F':147.068414,\
    'P':97.052764, 'S':87.032028, 'T':101.047679, 'U':150.95363, 'W':186.079313,\
    'Y':163.06332, 'V':99.068414, 'H2O':18.01056,'Proton':1.0072766}

mpModMass = {
    'Carbamidomethyl[C]': 57.021464,
    'Oxidation[M]': 15.994915
}  # 添加的修饰名称及质量
mstol = 20.0 


def calMassPepPlusMod(massList, modMass):
    pepMass = sum(massList) + mpMassTable['H2O']
    return  pepMass +  modMass + mpMassTable['Proton']


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
                pos -= 1             # 与python的index 从0开始对齐
                if pos == pep_a_len: #考虑修饰在C端的情况此时pos = peptide——A长度
                    pos -= 1         #超出了list的index，因此减一
                elif pos == -1:      # 考虑修饰在N端的情况 pos = -1，因此+1，算到第一个氨基酸是
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

    ####regular b,y ions ####
    pep_mass = a_long_mass + b_short_mass + mpMassTable['H2O']
    a_linker_mass = sum(a_mass_list) + mpMassTable['H2O'] + LinkerMass
    b_linker_mass = sum(b_mass_list) + mpMassTable['H2O'] + LinkerMass

    tmp_mass = mpMassTable['Proton']
    # b1+ of alpha
    b1Ion_pep_a = [0] * (pep_a_len - 1)  
    for i in range(pep_a_len - 1):
        if i == site_a:
            tmp_mass += b_linker_mass
        else:
            pass

        tmp_mass += a_mass_list[i]
        b1Ion_pep_a[i] = tmp_mass

    tmp_mass = mpMassTable['Proton'] + mpMassTable['H2O']
    y1Ion_pep_a = [0] * (pep_a_len - 1)  # y1+ of alpha
    for i in range(pep_a_len - 2, -1, -1):
        tmp_mass += a_mass_list[i + 1]
        if i + 1 == site_a:
            tmp_mass += b_linker_mass
        y1Ion_pep_a[i] = tmp_mass

    tmp_mass = mpMassTable['Proton']
    b1_b = [0] * (pep_b_len - 1)  # b1+ of beta
    for i in range(pep_b_len - 1):
        tmp_mass += b_mass_list[i]
        if i == site_b:
            tmp_mass += a_linker_mass
        b1_b[i] = tmp_mass

    tmp_mass = mpMassTable['Proton'] + mpMassTable['H2O']
    y1_b = [0] * (pep_b_len - 1)  # y1+ of beta
    for i in range(pep_b_len - 2, -1, -1):
        tmp_mass += b_mass_list[i + 1]
        if i + 1 == site_b:
            tmp_mass += a_linker_mass
        y1_b[i] = tmp_mass

    #######################Internal ions mass 1+  cleave at linker and backbone####################
    In_pepA_long = [0] * (pep_a_len - 1)
    for i in range(pep_a_len - 1):
        if i < site_a:
            In_pepA_long[i] = sum(
                a_mass_list[i + 1:]) + mpMassTable['H2O'] + longArmMass + 1.0078
        else:
            In_pepA_long[i] = sum(a_mass_list[:i + 1]) + longArmMass + 1.0078

    In_pepA_short = [0] * (pep_a_len - 1)
    for i in range(pep_a_len - 1):
        if i < site_a:
            In_pepA_short[i] = sum(
                a_mass_list[i + 1:]) + mpMassTable['H2O'] + shortArmMass + 1.0078
        else:
            In_pepA_short[i] = sum(a_mass_list[:i + 1]) + shortArmMass + 1.0078

    In_pepB_long = [0] * (pep_b_len - 1)
    for i in range(pep_b_len - 1):
        if i < site_b:
            In_pepB_long[i] = sum(
                b_mass_list[i + 1:]) + mpMassTable['H2O'] + longArmMass + 1.0078
        else:
            In_pepB_long[i] = sum(b_mass_list[:i + 1]) + longArmMass + 1.0078

    In_pepB_short = [0] * (pep_b_len - 1)
    for i in range(pep_b_len - 1):
        if i < site_b:
            In_pepB_short[i] = sum(
                b_mass_list[i + 1:]) + mpMassTable['H2O'] + shortArmMass + 1.0078
        else:
            In_pepB_short[i] = sum(b_mass_list[:i + 1]) + shortArmMass + 1.0078

    #######################Internal ions mass 1+  cleave at linker and backbone####################

    return [b1Ion_pep_a, y1Ion_pep_a, b1_b,
            y1_b], [a_long_mass, a_short_mass, b_long_mass, b_short_mass],\
                [In_pepA_long, In_pepA_short, In_pepB_long, In_pepB_short]


def getTargetSpec(spec_path):
    f = open(spec_path).readlines()
    scanSpecdic = {}
    i = 0
    while i < len(f):
        if f[i].strip() != "BEGIN IONS":
            print(i)
        else:
            scan = int(f[i+2].split("=")[1].strip())
            vOneSpec = []
            p = i+6
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

print(getTargetSpec(mgfPath))


def isMatched(mz, spec):
    for i in range(len(spec)):
        if 'END' in spec[i]:
            break
        exp_mz = float(spec[i].split(' ')[0])
        if abs((exp_mz - mz) / mz) * 1000000.0 <= mstol:
            return True
    return False


def generate_ion_mass_range(num):
    deta = num * mstol / 1000000
    return num - deta, num + deta

def match_report(mass_list, spec, max_c):
    ms2_dic = {}
    for i in range(len(spec) - 1):
        ms2_mz = float(spec[i][:-1].split(" ")[0])
        ms2_ins = float(spec[i][:-1].split(" ")[1])
        ms2_dic[ms2_mz] = ms2_ins
    rep_mask = [0] * 4
    match_dic = {}
    mz_list = sorted(list(ms2_dic.keys()))
    # print(mz_list)
    max_ins = max(list(ms2_dic.values()))
    for c in range(1, max_c + 1):
        mass_cur_list = [(m + 1.0078 * (c - 1)) / c for m in mass_list]

        repo_up_list = []
        repo_range_dic = {}
        for mass in mass_cur_list:
            low_mass, high_mass = generate_ion_mass_range(mass)
            repo_up_list.append(high_mass)
            repo_range_dic[mass] = [low_mass, high_mass]

            ft_idx = mass_cur_list.index(mass)
            ms2_idx = 0
            while ms2_idx < len(mz_list):
                if mz_list[ms2_idx] > high_mass:
                    break
                elif mz_list[ms2_idx] >= low_mass:
                    rep_mask[ft_idx] += 1
                    if ft_idx not in match_dic:
                        match_dic[ft_idx] = [[
                            c, mass_cur_list[ft_idx], mz_list[ms2_idx],
                            ms2_dic[mz_list[ms2_idx]] / max_ins
                        ]]
                    else:
                        match_dic[ft_idx].append([
                            c, mass_cur_list[ft_idx], mz_list[ms2_idx],
                            ms2_dic[mz_list[ms2_idx]] / max_ins
                        ])
                    break
                else:
                    ms2_idx += 1
    return match_dic


def summary_rep_macthDic(match_dic):
    pairBeta_list = []
    if 0 in match_dic and 1 in match_dic:
        chrg0_list = []
        chrg1_list = []
        ints0_list = []
        ints1_list = []
        for wd in match_dic[0]:
            chrg0_list.append(wd[0])
            ints0_list.append(wd[-1])
        for wd in match_dic[1]:
            chrg1_list.append(wd[0])
            ints1_list.append(wd[-1])
        xBetaPair_bool = False
        if len(chrg0_list) <= len(chrg1_list):
            for chg in chrg0_list:
                if chg in chrg1_list:
                    xBetaPair_bool = True
                    ave_ints = (ints0_list[chrg0_list.index(chg)] +
                                ints1_list[chrg1_list.index(chg)]) / 2
                    pairBeta_list.append([chg, ave_ints])
        else:
            for chg in chrg1_list:
                if chg in chrg0_list:
                    xBetaPair_bool = True
                    ave_ints = (float(ints0_list[chrg0_list.index(chg)]) +
                                float(ints1_list[chrg1_list.index(chg)])) / 2
                    pairBeta_list.append([chg, ave_ints])
    else:
        xBetaPair_bool = False

    pairAlpha_list = []
    if 2 in match_dic and 3 in match_dic:
        chrg2_list = []
        chrg3_list = []
        ints2_list = []
        ints3_list = []
        for wd in match_dic[2]:
            chrg2_list.append(wd[0])
            ints2_list.append(wd[-1])
        for wd in match_dic[3]:
            chrg3_list.append(wd[0])
            ints3_list.append(wd[-1])
        xAlphapair_bool = False
        if len(chrg2_list) <= len(chrg3_list):
            for chg in chrg2_list:
                if chg in chrg3_list:
                    xAlphapair_bool = True
                    ave_ints = (float(ints2_list[chrg2_list.index(chg)]) +
                                float(ints3_list[chrg3_list.index(chg)])) / 2
                    pairAlpha_list.append([chg, ave_ints])
        else:
            for chg in chrg3_list:
                if chg in chrg2_list:
                    xAlphapair_bool = True
                    ave_ints = (float(ints2_list[chrg2_list.index(chg)]) +
                                float(ints3_list[chrg3_list.index(chg)])) / 2
                    pairAlpha_list.append([chg, ave_ints])
    else:
        xAlphapair_bool = False

    detect_bool = False
    pair_num = [xBetaPair_bool, xAlphapair_bool].count(True)
    if pair_num > 0:
        detect_bool = True
    rep_ints_beta = "na"
    if pairBeta_list:
        pairBeta_ints_list = []
        for pr in pairBeta_list:
            pairBeta_ints_list.append(pr[-1])
        rep_ints_beta = str(max(pairBeta_ints_list))

    rep_ints_alpha = "na"
    if pairAlpha_list:
        pairalpha_ints_list = []
        for pr in pairAlpha_list:
            pairalpha_ints_list.append(pr[-1])
        rep_ints_alpha = str(max(pairalpha_ints_list))

    return detect_bool, pair_num, rep_ints_beta, rep_ints_alpha


def cal1pep(common_mask, xlink_mask, b1, y1, site, max_c, length, spec):
    for c in range(1, max_c + 1):
        b_cur = []
        y_cur = []
        if c == 1:
            b_cur = copy.deepcopy(b1)  # 1+ b ion
            y_cur = copy.deepcopy(y1)  # 1+ y ion
        else:
            b_cur = [(m + mpMassTable['Proton'] * (c - 1)) / float(c)
                     for m in b1]
            y_cur = [(m + mpMassTable['Proton'] * (c - 1)) / float(c)
                     for m in y1]

        for i in range(length - 1):
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


def cal1pepReporter(a_long_mass, a_short_mass, b_long_mass, b_short_mass,
                    max_c, spec):
    reP_mask = zeros([4, 1])
    alm, asm, blm, bsm = a_long_mass, a_short_mass, b_long_mass, b_short_mass

    for c in range(1, max_c + 1):
        alm = (alm + mpMassTable['Proton'] * (c - 1)) / float(c)
        asm = (asm + mpMassTable['Proton'] * (c - 1)) / float(c)
        blm = (blm + mpMassTable['Proton'] * (c - 1)) / float(c)
        bsm = (bsm + mpMassTable['Proton'] * (c - 1)) / float(c)

        if isMatched(alm, spec):
            fmatched.write('charge=%d\tmz=%f\ttype=a_long_mass\n' % (c, alm))
            reP_mask[0] += 1

        if isMatched(asm, spec):
            fmatched.write('charge=%d\tmz=%f\ttype=a_short_mass\n' % (c, asm))
            reP_mask[1] += 1

        if isMatched(blm, spec):
            fmatched.write('charge=%d\tmz=%f\ttype=b_long_mass\n' % (c, blm))
            reP_mask[2] += 1

        if isMatched(bsm, spec):
            fmatched.write('charge=%d\tmz=%f\ttype=b_short_mass\n' % (c, bsm))
            reP_mask[3] += 1
    return reP_mask


def cal1pepIntXP(In_pepA_long, In_pepA_short, In_pepB_long, In_pepB_short,spec):
    max_c = 2
    inPAXL = copy.deepcopy(In_pepA_long)
    inPAXS = copy.deepcopy(In_pepA_short)
    inPBXL = copy.deepcopy(In_pepB_long)
    inPBXS = copy.deepcopy(In_pepB_short)
    inPAXL_mask = [0] * len(inPAXL)
    inPAXS_mask = [0] * len(inPAXS)
    inPBXL_mask = [0] * len(inPBXL)
    inPBXS_mask = [0] * len(inPBXS)

    for c in range(1, max_c + 1):
        inPAXL_cur = [(m + 1.0078 * (c - 1)) / c for m in inPAXL]
        inPAXS_cur = [(m + 1.0078 * (c - 1)) / c for m in inPAXS]
        inPBXL_cur = [(m + 1.0078 * (c - 1)) / c for m in inPBXL]
        inPBXS_cur = [(m + 1.0078 * (c - 1)) / c for m in inPBXS]

        for i in range(len(inPAXL_cur)):
            if isMatched(inPAXL_cur[i], spec):
                inPAXL_mask[i] += 1
                fmatched.write(
                    'site=(%d,AL)\tcharge=%d\tmz=%f\ttype=In_pepA_long\n' %
                    (i, c, inPAXL_cur[i]))
            if isMatched(inPAXS_cur[i], spec):
                inPAXS_mask[i] += 1
                fmatched.write(
                    'site=(%d,AS)\tcharge=%d\tmz=%f\ttype=In_pepA_short\n' %
                    (i, c, inPAXS_cur[i]))

        for i in range(len(inPBXL_cur)):
            if isMatched(inPBXL_cur[i], spec):
                inPBXL_mask[i] += 1
                fmatched.write(
                    'site=(%d,BL)\tcharge=%d\tmz=%f\tIn_pepB_long\n' %
                    (i, c, inPBXL_cur[i]))
            if isMatched(inPBXS_cur[i], spec):
                inPBXS_mask[i] += 1
                fmatched.write(
                    'site=(%d,BS)\tcharge=%d\tmz=%f\tIn_pepB_short\n' %
                    (i, c, inPBXS_cur[i]))

    total_unmat = inPAXL_mask.count(0) + inPAXS_mask.count(
        0) + inPBXL_mask.count(0) + inPBXS_mask.count(0)
    intXP_ratio = 1 - total_unmat / float(len(inPAXL) + len(inPBXL)) / 2
    return intXP_ratio  # , [inPAXL_mask, inPAXS_mask, inPBXL_mask, inPBXS_mask]


def calIonRatio(pep, spec, charge, modCell):

    a, b = pep.split('-')

    site_a = int(''.join([i for i in a if i.isdigit()])) - 1
    site_b = int(''.join([i for i in b if i.isdigit()])) - 1

    a = ''.join([i for i in a if i.isalpha()])
    b = ''.join([i for i in b if i.isalpha()])

    pep_a_len = len(a)
    pep_b_len = len(b)

    a_common_mask = [0] * (pep_a_len - 1)
    a_xlink_mask = [0] * (pep_a_len - 1)

    b_common_mask = [0] * (pep_b_len - 1)
    b_xlink_mask = [0] * (pep_b_len - 1)

    by_ions, reporter_ions, intXPions = getTheroIons(pep, modCell)

    b1_a, y1_a, b1_b, y1_b = by_ions
    a_long_mass, a_short_mass, b_long_mass, b_short_mass = reporter_ions

    In_pepA_long, In_pepA_short, In_pepB_long, In_pepB_short = intXPions

    ftheo_ions.write('title=%s\tcharge=%s\n' % (pep, charge))
    
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

    max_c = min(3, charge)

    fmatched.write('title=%s\n' % (pep))
    fmatched.write('---------alpha---------------\n')
    cal1pep(a_common_mask, a_xlink_mask, b1_a, y1_a, site_a, max_c, pep_a_len,
            spec)
    fmatched.write('common ion matched: ' + str(a_common_mask) + '\n')
    fmatched.write('xlink ion matched: ' + str(a_xlink_mask) + '\n')

    fmatched.write('---------beta---------------\n')
    cal1pep(b_common_mask, b_xlink_mask, b1_b, y1_b, site_b, max_c, pep_b_len,
            spec)
    fmatched.write('common ion mached: ' + str(b_common_mask) + '\n')
    fmatched.write('xlink ion matched: ' + str(b_xlink_mask) + '\n')

    fmatched.write('---------reporter ions---------------\n')
    cal1pepReporter(a_long_mass, a_short_mass, b_long_mass, b_short_mass,
                    max_c, spec)
    report_mactch_dic = match_report(
        [b_short_mass, b_long_mass, a_short_mass, a_long_mass], spec, max_c)
    detect_bool, pair_num, rep_ints_beta, rep_ints_alpha = summary_rep_macthDic(
        report_mactch_dic)

    fmatched.write('---------internal PX---------------\n')
    intPXratio = cal1pepIntXP(In_pepA_long, In_pepA_short, In_pepB_long,
                           In_pepB_short, spec)

    sumCommAlpha = sum(a_common_mask)
    sumXlinkAlpha = sum(a_xlink_mask)
    sumCommBeta = sum(b_common_mask)
    sumXlinkBeta = sum(b_xlink_mask)

    regular_by_ion_ratio = (sumCommAlpha + sumCommBeta) / float(pep_b_len + pep_a_len - 2)
    xlink_by_ion_ratio = (sumXlinkAlpha + sumXlinkBeta) / float(pep_b_len + pep_a_len - 2)
    commByAlpha = sumCommAlpha/(pep_a_len -1)
    commByBeta = sumCommBeta/(pep_b_len - 1)
    
    return regular_by_ion_ratio, xlink_by_ion_ratio, detect_bool, pair_num, rep_ints_beta, rep_ints_alpha, intPXratio, commByAlpha, commByBeta


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
            #print(scan)
            mz = float(f[i].split("\t")[-1])
            charge = f[i+9].split("\t")[1]
            nce = f[i+6].split("\t")[2].split(" ")[7].split("@")[1][3:5]
            #print(charge)
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


scanNCEdic = linkscanNCEmzPepMod()

repDic = {}
scanSpecdic = getTargetSpec(mgfPath)
for scan in scanNCEdic:
    ttl = scanNCEdic[scan]
    #print(scan, ttl)
    pep = ttl[3]
    modCell = ttl[4]
    charge = int(ttl[1])
    nce = ttl[2]
    
    spec = scanSpecdic[scan]

    extendInfo = calIonRatio(pep, spec, charge, modCell)
    ttl.extend(extendInfo)
    ttl.insert(0, scan)
    repDic[scan] = ttl

fmatched.close()
ftheo_ions.close()

b = open("reportFile_commBY.csv", 'w')
b.write("scan, mz, charge, nce, pep, mods, regular_by_ion_ratio,\
     xlink_by_ion_ratio, detect_bool, pair_num, rep_ints_beta,\
         rep_ints_alpha, intPXratio, commBYalpha, commBYbeta" + "\n")
#nceList = sorted(list(repDic.keys()))
#for nce in nceList:
#scanList = 
for scan in repDic:
    wlist = [str(ele) for ele in repDic[scan]]
    b.write(",".join(wlist)+"\n")
b.close()

main()