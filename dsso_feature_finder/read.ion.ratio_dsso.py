# -*- coding: utf-8 -*-
"""
Created on Mon May 07 09:43:52 2018

@author: zhenlin
"""
import os
import copy

# mgf 所在文件夹路径
mgf_dir = r'E:\Collaboration\NIBS\caoyong\CXMS_opt\example\mgf'
# pLink2 鉴定到的交联结果文件路径，谱图层次
id_path = r'E:\Collaboration\NIBS\caoyong\CXMS_opt\example\pLink2_id_results\uniprot-ecoli-20171023_2017.12.22.filtered_cross-linked_spectra.csv'
# 生成的匹配信息文件路径
matched_info = r'E:\Collaboration\NIBS\caoyong\CXMS_opt\example\output\matched_info.txt'
# 生成的Ion Ratio信息文件路径
ion_ratio_info = r'E:\Collaboration\NIBS\caoyong\CXMS_opt\example\output\ion_ratio_info.csv'
# 输出所有理论离子
theoretical_ions = r'E:\Collaboration\NIBS\caoyong\CXMS_opt\example\output\theoretical_ions.txt'

LinkerMass = 158.004  # 交联剂质量
LongMass = 85.9826
ShortMass = 54.0106

mpModMass = {
    'Carbamidomethyl[C]': 57.021464,
    'Oxidation[M]': 15.994915
}  # 添加的修饰名称及质量
mstol = 20.0  # 碎片离子匹配误差，ppm为单位

#########################  请勿修改以下内容  ########################

mpMassTable={'A':71.037114,'R':156.101111,'N':114.042927,'D':115.026943,'C':103.009185, \
'E':129.042593,'Q':128.058578,'G':57.021464,'H':137.058912,'I':113.084064, \
'L':113.084064,'K':128.094963,'M':131.040485,'F':147.068414,'P':97.052764, \
'S':87.032028,'T':101.047679,'U':150.95363,'W':186.079313,'Y':163.06332,'V':99.068414, \
'H2O':18.01056,'Proton':1.0072766}

fmatched = open(matched_info, 'w')
ftheo_ions = open(theoretical_ions, 'w')


# 读取所有psm信息
def getAllPSM(path):
    fid = open(path)
    lines = fid.readlines()
    fid.close()

    mpPSM = {}
    for i in range(1, len(lines)):
        lines[i] = lines[i].strip()
        if len(lines[i]) < 1:
            continue
        contents = lines[i].split(',')

        order = int(contents[0])
        title = contents[1]
        charge = int(contents[2])
        peptide = contents[4]
        mods = contents[8]

        psm = [
            order, title, charge, peptide, mods, 0, 0, 0, 0, 0, 0
        ]  # a_len, b_len, alpha common b/y mic, alpha xlink b/y mic, beta common b/y mic, beta xlink b/y mic

        raw = title.split('.')[0]
        if raw in mpPSM:
            mpPSM[raw].append(psm)
        else:
            mpPSM[raw] = [psm]
    return mpPSM


# 针对一个psm，计算理论碎片离子，都是1价的
def getIons(psm):
    a, b = psm[3].split('-')

    a_site = int(''.join([i for i in a if i.isdigit()])) - 1
    b_site = int(''.join([i for i in b if i.isdigit()])) - 1

    a = ''.join([i for i in a if i.isalpha()])
    b = ''.join([i for i in b if i.isalpha()])

    a_len = len(a)
    b_len = len(b)

    a_mod_mass = [0] * a_len # a modifiaction_mass_list
    b_mod_mass = [0] * b_len # b modi_mass_list

    if psm[4] != 'null':
        mods = psm[4].split(';')
        for m in mods:
            name, pos = m.split('(')
            pos = int(''.join([i for i in pos if i.isdigit()]))
            if pos <= a_len + 1:
                pos -= 1
                if pos == a_len:
                    pos -= 1
                a_mod_mass[pos] += mpModMass[name]
            else:
                pos = pos - a_len - 4
                if pos == -1:
                    pos = 0
                b_mod_mass[pos] += mpModMass[name]
    else:
        pass

    a_mass = [0] * a_len  # 残基质量
    b_mass = [0] * b_len

    for i in range(a_len):
        a_mass[i] = mpMassTable[a[i]] + a_mod_mass[i]
    for i in range(b_len):
        b_mass[i] = mpMassTable[b[i]] + b_mod_mass[i]

# ####################reporter ions###########################################
    a_long_mass = sum(
        a_mass) + mpMassTable['H2O'] + LongMass + mpMassTable['Proton']
    a_short_mass = sum(
        a_mass) + mpMassTable['H2O'] + ShortMass + mpMassTable['Proton']
    b_long_mass = sum(
        b_mass) + mpMassTable['H2O'] + LongMass + mpMassTable['Proton']
    b_short_mass = sum(
        b_mass) + mpMassTable['H2O'] + ShortMass + mpMassTable['Proton']
# ####################reporter ions###########################################

    a_linker_mass = sum(a_mass) + mpMassTable['H2O'] + LinkerMass
    b_linker_mass = sum(b_mass) + mpMassTable['H2O'] + LinkerMass

    tmp_mass = mpMassTable['Proton']
    b1_a = [0] * (a_len - 1)  # b1+ of alpha
    for i in range(a_len - 1):
        tmp_mass += a_mass[i]
        if i == a_site:
            tmp_mass += b_linker_mass
        b1_a[i] = tmp_mass

    tmp_mass = mpMassTable['Proton'] + mpMassTable['H2O']
    y1_a = [0] * (a_len - 1)  # y1+ of alpha
    for i in range(a_len - 2, -1, -1):
        tmp_mass += a_mass[i + 1]
        if i + 1 == a_site:
            tmp_mass += b_linker_mass
        y1_a[i] = tmp_mass

    tmp_mass = mpMassTable['Proton']
    b1_b = [0] * (b_len - 1)  # b1+ of beta
    for i in range(b_len - 1):
        tmp_mass += b_mass[i]
        if i == b_site:
            tmp_mass += a_linker_mass
        b1_b[i] = tmp_mass

    tmp_mass = mpMassTable['Proton'] + mpMassTable['H2O']
    y1_b = [0] * (b_len - 1)  # y1+ of beta
    for i in range(b_len - 2, -1, -1):
        tmp_mass += b_mass[i + 1]
        if i + 1 == b_site:
            tmp_mass += a_linker_mass
        y1_b[i] = tmp_mass

#####################internal ions###########################################
    a_mass[a_site] += LongMass
    a_long_by = [[0] * a_len] * a_len

    for i in range(0, a_len):
        tmp_mass = mpMassTable['Proton']
        for j in range(i, a_len):
            tmp_mass += a_mass[j]
            a_long_by[i][j] = tmp_mass

    a_mass[a_site] = a_mass[a_site] - LongMass + ShortMass
    a_short_by = [[0] * a_len] * a_len

    for i in range(0, a_len):
        tmp_mass = mpMassTable['Proton']
        for j in range(i, a_len):
            tmp_mass += a_mass[j]
            a_short_by[i][j] = tmp_mass

    b_mass[b_site] += LongMass
    b_long_by = [[0] * b_len] * b_len

    for i in range(0, b_len):
        tmp_mass = mpMassTable['Proton']
        for j in range(i, b_len):
            tmp_mass += b_mass[j]
            b_long_by[i][j] = tmp_mass

    b_mass[b_site] = b_mass[b_site] - LongMass + ShortMass
    b_short_by = [[0] * b_len] * b_len

    for i in range(0, b_len):
        tmp_mass = mpMassTable['Proton']
        for j in range(i, b_len):
            tmp_mass += b_mass[j]
            b_short_by[i][j] = tmp_mass

#####################internal ions###########################################

    return [b1_a, y1_a, b1_b, y1_b], [
        a_long_mass, a_short_mass, b_long_mass, b_short_mass
    ], [a_long_by, a_short_by, b_long_by, b_short_by]


# 从mgf文件中提取需要用到的谱图
def getTargetSpec(PSMs, spec_path):
    titles = set()
    for psm in PSMs:
        titles.add(psm[1])

    mpTitleSpec = {}

    fIn = open(spec_path)

    vOneSpec = []
    for line in fIn:
        if 'BEGIN' in line:
            if len(vOneSpec) > 0 and vOneSpec[1] in titles:
                mpTitleSpec[vOneSpec[1]] = copy.deepcopy(vOneSpec)
            vOneSpec = []
            vOneSpec.append([line])
        else:
            vOneSpec[0].append(line)
        if 'TITLE' in line:
            title = line.strip().split('=')[-1].strip()
            vOneSpec.append(title)
    if len(vOneSpec) > 0 and vOneSpec[1] in titles:
        mpTitleSpec[vOneSpec[1]] = copy.deepcopy(vOneSpec)
    return mpTitleSpec


def isMatched(mz, spec):
    for i in range(5, len(spec[0])):
        if 'END' in spec[0][i]:
            break
        exp_mz = float(spec[0][i].split(' ')[0])
        if abs((exp_mz - mz) / mz) * 1000000.0 <= mstol:
            return True
    return False


def cal1pep(common_mask, xlink_mask, b1, y1, site, max_c, length, spec):

    for c in range(1, max_c + 1):
        b_cur = []
        y_cur = []
        if c == 1:
            b_cur = copy.deepcopy(b1)
            y_cur = copy.deepcopy(y1)
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

    alm, asm, blm, bsm = a_long_mass, a_short_mass, b_long_mass, b_short_mass

    for c in range(1, max_c + 1):
        alm = (alm + mpMassTable['Proton'] * (c - 1)) / float(c)
        asm = (asm + mpMassTable['Proton'] * (c - 1)) / float(c)
        blm = (blm + mpMassTable['Proton'] * (c - 1)) / float(c)
        bsm = (bsm + mpMassTable['Proton'] * (c - 1)) / float(c)

        if isMatched(alm, spec):
            fmatched.write('charge=%d\tmz=%f\ttype=a_long_mass\n' % (c, alm))

        if isMatched(asm, spec):
            fmatched.write('charge=%d\tmz=%f\ttype=a_short_mass\n' % (c, asm))

        if isMatched(blm, spec):
            fmatched.write('charge=%d\tmz=%f\ttype=b_long_mass\n' % (c, blm))

        if isMatched(bsm, spec):
            fmatched.write('charge=%d\tmz=%f\ttype=b_short_mass\n' % (c, bsm))


def cal1pepInternal(a_long_by, a_short_by, b_long_by, b_short_by, max_c, spec):
    alby = copy.deepcopy(a_long_by)
    asby = copy.deepcopy(a_short_by)
    blby = copy.deepcopy(b_long_by)
    bsby = copy.deepcopy(b_short_by)

    for c in range(1, max_c + 1):

        for i in range(len(alby)):
            for j in range(len(alby[i])):
                mz = (alby[i][j] + mpMassTable['Proton'] * (c - 1)) / float(c)
                if isMatched(mz, spec):
                    fmatched.write(
                        'site=(%d,%d)\tcharge=%d\tmz=%f\ttype=a_long_by\n' %
                        (i, j, c, mz))

        for i in range(len(asby)):
            for j in range(len(asby[i])):
                mz = (asby[i][j] + mpMassTable['Proton'] * (c - 1)) / float(c)
                if isMatched(mz, spec):
                    fmatched.write(
                        'site=(%d,%d)\tcharge=%d\tmz=%f\ttype=a_short_by\n' %
                        (i, j, c, mz))

        for i in range(len(blby)):
            for j in range(len(blby[i])):
                mz = (blby[i][j] + mpMassTable['Proton'] * (c - 1)) / float(c)
                if isMatched(mz, spec):
                    fmatched.write(
                        'site=(%d,%d)\tcharge=%d\tmz=%f\ttype=b_long_by\n' %
                        (i, j, c, mz))

        for i in range(len(bsby)):
            for j in range(len(bsby[i])):
                mz = (bsby[i][j] + mpMassTable['Proton'] * (c - 1)) / float(c)
                if isMatched(mz, spec):
                    fmatched.write(
                        'site=(%d,%d)\tcharge=%d\tmz=%f\ttype=b_short_by\n' %
                        (i, j, c, mz))


def calIonRatio(psm, spec):

    a, b = psm[3].split('-')

    a_site = int(''.join([i for i in a if i.isdigit()])) - 1
    b_site = int(''.join([i for i in b if i.isdigit()])) - 1

    a = ''.join([i for i in a if i.isalpha()])
    b = ''.join([i for i in b if i.isalpha()])

    a_len = len(a)
    b_len = len(b)
    psm[5] = a_len
    psm[6] = b_len

    a_common_mask = [0] * (a_len - 1)
    a_xlink_mask = [0] * (a_len - 1)

    b_common_mask = [0] * (b_len - 1)
    b_xlink_mask = [0] * (b_len - 1)

    by_ions, reporter_ions, internal_ions = getIons(psm)

    b1_a, y1_a, b1_b, y1_b = by_ions
    a_long_mass, a_short_mass, b_long_mass, b_short_mass = reporter_ions
    a_long_by, a_short_by, b_long_by, b_short_by = internal_ions

    ftheo_ions.write('title=%s\n' % (psm[1]))
    ftheo_ions.write(psm + '\n')

    ftheo_ions.write('--------alpha backbone ions (alpha b1+)---------\n')
    ftheo_ions.write(b1_a)

    ftheo_ions.write('--------alpha backbone ions (alpha y1+)---------\n')
    ftheo_ions.write(y1_a)

    ftheo_ions.write('--------beta backbone ions (beta b1+)---------\n')
    ftheo_ions.write(b1_b)

    ftheo_ions.write('--------beta backbone ions (beta y1+)---------\n')
    ftheo_ions.write(y1_b)

    ftheo_ions.write('--------reporter ions (alpha long 1+)---------\n')
    ftheo_ions.write(a_long_mass)

    ftheo_ions.write('--------reporter ions (alpha short 1+)---------\n')
    ftheo_ions.write(a_short_mass)

    ftheo_ions.write('--------reporter ions (beta long 1+)---------\n')
    ftheo_ions.write(b_long_mass)

    ftheo_ions.write('--------reporter ions (beta short 1+)---------\n')
    ftheo_ions.write(b_short_mass)

    ftheo_ions.write('--------internal ions (alpha long by 1+)---------\n')
    ftheo_ions.write(a_long_by)

    ftheo_ions.write('--------internal ions (alpha short by 1+)---------\n')
    ftheo_ions.write(a_short_by)

    ftheo_ions.write('--------internal ions (beta long by 1+)---------\n')
    ftheo_ions.write(b_long_by)

    ftheo_ions.write('--------internal ions (beta short by 1+)---------\n')
    ftheo_ions.write(b_short_by)

    max_c = min(3, psm[2])

    fmatched.write('title=%s\n' % (psm[1]))
    fmatched.write('---------alpha---------------\n')
    cal1pep(a_common_mask, a_xlink_mask, b1_a, y1_a, a_site, max_c, a_len,
            spec)
    fmatched.write('common ion matched: ' + str(a_common_mask) + '\n')
    fmatched.write('xlink ion matched: ' + str(a_xlink_mask) + '\n')

    fmatched.write('---------beta---------------\n')
    cal1pep(b_common_mask, b_xlink_mask, b1_b, y1_b, b_site, max_c, b_len,
            spec)
    fmatched.write('common ion mached: ' + str(b_common_mask) + '\n')
    fmatched.write('xlink ion matched: ' + str(b_xlink_mask) + '\n')

    fmatched.write('---------reporter ions---------------\n')
    cal1pepReporter(a_long_mass, a_short_mass, b_long_mass, b_short_mass,
                    max_c, spec)

    fmatched.write('---------internal ions---------------\n')
    cal1pepInternal(a_long_by, a_short_by, b_long_by, b_short_by, max_c, spec)

    psm[7] = sum(a_common_mask)
    psm[8] = sum(a_xlink_mask)
    psm[9] = sum(b_common_mask)
    psm[10] = sum(b_xlink_mask)


mpPSM = getAllPSM(id_path)

ans = []
for raw, PSMs in mpPSM.items():
    mpTitleSpec = getTargetSpec(PSMs, os.path.join(mgf_dir,
                                                   raw + '_HCDFT.mgf'))
    for psm in PSMs:
        calIonRatio(psm, mpTitleSpec[psm[1]])
        ans.append(psm)

fmatched.close()

ans = sorted(ans, key=lambda x: x[0])

fionratio = open(ion_ratio_info, 'w')
fionratio.write('Order,Title,Charge,Peptide,Modification,alpha length, beta length, regular b/y matched count of alpha, xlink b/y matched count of alpha,' \
'regular b/y matched count of beta, xlink b/y matched count of beta, ion ratio of regular b/y, ion ratio of xlink b/y\n')
for p in ans:
    regular_by_ion_ratio = (p[7] + p[9]) / float(p[5] + p[6] - 2)
    xlink_by_ion_ratio = (p[8] + p[10]) / float(p[5] + p[6] - 2)
    fionratio.write('%d,%s,%d,%s,%s,%d,%d,%d,%d,%d,%d,%f,%f\n' %
                    (p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8],
                     p[9], p[10], regular_by_ion_ratio, xlink_by_ion_ratio))
fionratio.close()