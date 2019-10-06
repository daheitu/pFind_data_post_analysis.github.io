# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:52:57 2019

@author: Zhenlin
"""
import os
from scipy.stats import norm

path = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.10.02.18.35.44_Lacto_DSSO_R3"

pLink2_id_folder = path
inclusion_list = os.path.join(path, 'inclusion_list.csv')
inclusion_list_detail = os.path.join(path, 'inclusion_list_detail.csv')

######## do not touch files below  #####################
pq = r'%s\pQuant\pQuant.spectra' % pLink2_id_folder

MP = 1.0072766

specid = ''
param = ''
mgfs = []
ms1s = []
raws = []

reports = os.listdir(r'%s\reports' % pLink2_id_folder)
for c in reports:
    if 'filtered_cross-linked_spectra.csv' in c:
        specid = r'%s\reports\%s' % (pLink2_id_folder, c)

roots = os.listdir(pLink2_id_folder)
for c in roots:
    if '.plink' in c:
        param = r'%s\%s' % (pLink2_id_folder, c)

fp = open(param)
lines = fp.readlines()
fp.close()

for i in range(len(lines)):
    if '[spectrum]' in lines[i]:
        for j in range(i, len(lines)):
            if 'spec_path' in lines[j]:
                p = lines[j].split('=')[1].strip()
                p = r'%s' % p

                fd = os.path.dirname(p)
                pf2 = os.path.basename(p)
                rawname = pf2[:pf2.index('HCDFT') - 1]

                raws.append(rawname)
                mgfs.append(r'%s\%s_HCDFT.mgf' % (fd, rawname))
                ms1s.append(r'%s\%s.ms1' % (fd, rawname))


# 读取mgf文件，返回{title:spectrum}
def ExtractMGFTitle3(strPathIn):
    fIn = open(strPathIn)
    vOneSpec = []
    ans = {}
    for line in fIn:
        if 'BEGIN' in line:
            if len(vOneSpec) > 0:
                ans[vOneSpec[1]] = vOneSpec
            vOneSpec = []
            vOneSpec.append([line])
        elif not line[0].isdigit():
            vOneSpec[0].append(line)
        if 'TITLE' in line:
            title = line.strip().split('=')[-1].strip()
            vOneSpec.append(title)
    if len(vOneSpec) > 0:
        ans[vOneSpec[1]] = vOneSpec
    fIn.close()
    return ans


# 读取ms1，输出dict, key=scan, value=rt
def readMS1(path):  #z
    fin = open(path)
    lines = fin.readlines()
    fin.close()
    mpScanRT = {}
    scan = ''
    rt = ''
    for line in lines:
        if line[0] == 'S':
            scan = line.strip().split('\t')[1]
            scan = int(scan)
        elif 'RetTime' in line:
            rt = line.strip().split('\t')[-1]
            rt = float(rt)
            mpScanRT[scan] = rt
    return mpScanRT


specs = {}
for mgf in mgfs:
    s = ExtractMGFTitle3(mgf)
    specs = {**specs, **s}

ms1specs = {}
for i in range(len(ms1s)):
    pms1 = ms1s[i]
    rawname = raws[i]
    s = readMS1(pms1)
    ms1specs[rawname] = s


# 读取pQuant.spectra
# 返回dict，key=title, value色谱曲线的信息=[[ms1 scan list],[mono intensity list], [rt-min,rt-max]]
def ReadpQuant(path):
    fin = open(path)
    lines = fin.readlines()
    fin.close()

    mpTitleChrom = {}

    n = len(lines)
    i = 0
    while i < n:
        title = lines[i].strip().split('\t')[-1]
        j = i + 1
        while j < n and lines[j][0] != '@':
            j += 1
        scans = []
        mono_int = []
        rt = []
        for k in range(i + 1, j):
            if 'I,E,0,04' in lines[k]:
                scans = lines[k].strip().split('\t')[1:]
                scans = [int(x) for x in scans]
            elif 'I,E,0,05' in lines[k]:
                mono_int = lines[k + 1].strip().split('\t')
                mono_int = [float(x) for x in mono_int]
            elif 'I,E,0,08' in lines[k]:
                contents = lines[k].strip().split('\t')
                st_idx = int(contents[1])
                end_idx = int(contents[2])
                st_rt = float(contents[-2])
                end_rt = float(contents[-1])
                rt = [st_rt, end_rt]
                scans = scans[st_idx:end_idx + 1]
                mono_int = mono_int[st_idx:end_idx + 1]
        mpTitleChrom[title] = [scans, mono_int, rt]
        i = j
    return mpTitleChrom


mpTitleChrom = ReadpQuant(pq)

fin = open(specid)
lines = fin.readlines()
fin.close()

ans = {}

for i in range(1, len(lines)):
    contents = lines[i].split(',')
    title = contents[1].strip()

    charge = int(contents[2].strip())
    pep = contents[4].strip()
    mod = contents[8].strip()

    mz = (float(contents[7]) + (charge - 1) * MP) / charge

    if title not in specs:
        continue

    ms2_rt = specs[title][0][3]
    ms2_rt = float(ms2_rt.split('=')[1])
    svmscore = float(contents[10])

    evalue = float(contents[9])

    ms1_scans, ms1_mono_int, ms1_rt = mpTitleChrom[title]
    ms1_rtmin, ms1_rtmax = ms1_rt

    ms1_rts = []
    rawname = title.split('.')[0]
    for ms1scan in ms1_scans:
        ms1_rts.append(ms1specs[rawname][ms1scan])

    pp = '%s|%s|%d' % (pep, mod, charge)
    if pp in ans:
        ans[pp].append([
            mz, ms2_rt, evalue, svmscore, title, ms1_rtmin, ms1_rtmax,
            ms1_scans, ms1_mono_int, ms1_rts
        ])
    else:
        ans[pp] = [[
            mz, ms2_rt, evalue, svmscore, title, ms1_rtmin, ms1_rtmax,
            ms1_scans, ms1_mono_int, ms1_rts
        ]]

fout = open(inclusion_list, 'w')
fout.write(
    'pep_pair,charge,mods,theo_m/z,id_rt_min,id_rt_max,spec_cnt,evalue_best,svm_best,svm_worst,pq_rt_min,pq_rt_max,fit_chrom_mean,fit_chrom_std,mean-3std,mean+3std\n'
)

fout2 = open(inclusion_list_detail, 'w')
fout2.write(
    'pep_pair,charge,mods,theo_m/z,id_rt_min,id_rt_max,spec_cnt,evalue_best,svm_best,svm_worst,pq_rt_min,pq_rt_max,fit_chrom_mean,fit_chrom_std,mean-3std,mean+3std\n'
)
fout2.write('ms1 scans\n')
fout2.write('ms1_rt\n')
fout2.write('ms1_mono_intensity\n')

for k, v in ans.items():
    pep, mod, charge = k.split('|')
    theo_mz = v[0][0]
    id_rtmin = min(v, key=lambda x: x[1])[1]
    id_rtmax = max(v, key=lambda x: x[1])[1]

    evalue_best = min(v, key=lambda x: x[2])[2]

    svm_best = min(v, key=lambda x: x[3])[3]
    svm_worst = max(v, key=lambda x: x[3])[3]

    specn = len(v)

    ms1_rtmin = min(v, key=lambda x: x[5])[5]
    ms1_rtmax = max(v, key=lambda x: x[6])[6]

    ms1_scans = []
    ms1_mono_int = []
    ms1_rt = []

    for feature in v:
        s, m, r = feature[-3:]
        ms1_scans.extend(s)
        ms1_mono_int.extend(m)
        ms1_rt.extend(r)

    # 把(value,count)转换为一维的序列，以便拟合1D正态分布
    max_int = max(ms1_mono_int)
    freq = [1000 * i / max_int for i in ms1_mono_int]  #缩小到[0,1000]的范围
    freq = [int(i + 1) for i in freq]  # +1平滑
    data = []
    for i in range(len(ms1_rt)):
        for n in range(freq[i]):
            data.append(ms1_rt[i])
    #######################
    mean, std = norm.fit(data)

    fout.write('%s,%s,%s,%f,%f,%f,%d,%E,%E,%E,%f,%f,%f,%f,%f,%f\n' %
               (pep, charge, mod, theo_mz, id_rtmin, id_rtmax, specn,
                evalue_best, svm_best, svm_worst, ms1_rtmin, ms1_rtmax, mean,
                std, mean - 3 * std, mean + 3 * std))

    fout2.write('%s,%s,%s,%f,%f,%f,%d,%E,%E,%E,%f,%f,%f,%f,%f,%f\n' %
                (pep, charge, mod, theo_mz, id_rtmin, id_rtmax, specn,
                 evalue_best, svm_best, svm_worst, ms1_rtmin, ms1_rtmax, mean,
                 std, mean - 3 * std, mean + 3 * std))
    fout2.write(','.join([str(x) for x in ms1_scans]) + '\n')
    fout2.write(','.join([str(x) for x in ms1_rt]) + '\n')
    fout2.write(','.join([str(x) for x in ms1_mono_int]) + '\n')

    #break

fout.close()
fout2.close()