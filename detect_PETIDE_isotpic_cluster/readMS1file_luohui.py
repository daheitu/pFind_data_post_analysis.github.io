import os
import random
from detect_isotopics_20220405 import detect_isotopic_cluster
import gc

gc.enable()

os.chdir(r'Z:\MS_DATA\20220827')  # rawconvert到处ms1所在路径
tgt_ion_file = r"G:\test\cy_list.csv"   # 离子的列表

def generate_mass_range(num, delta_ppm):
    delta = num * delta_ppm / 1000000
    return num - delta, num + delta


# 寻找符合质量约束的离子
def findOneMZ1(mz, ms1MzList,  tol_ppm=10):
    k = 0
    ms1MzList.sort()
    low_ms, up_ms = generate_mass_range(mz, tol_ppm)
    while k < len(ms1MzList):
        if float(ms1MzList[k]) > up_ms:
            return False, -1
        elif float(ms1MzList[k]) >= low_ms:
            return True, ms1MzList[k]
        else:
            k += 1
    if k == len(ms1MzList):
        return False, -1


# 寻找强度最高的肽段feature
def detectPepFeat(ms1_info_dic, intensCutoff):
    finalPepFeatDic = {}
    for scan in ms1_info_dic:
        print(scan)
        currentMzIntens_list = ms1_info_dic[scan][1]
        currChgedIsoFeatrues = detect_isotopic_cluster(currentMzIntens_list)
        currChgedIsoDic = {}
        for fea_info in currChgedIsoFeatrues:
            mz,charge = fea_info[:2]
            if charge == 1 and mz < 500:
                continue
            else:
                currChgedIsoDic[(mz, charge)] = list(fea_info[2:])
                currChgedIsoDic[(mz, charge)].append(scan)
        
        for mz,charge in currChgedIsoDic:
            intens = currChgedIsoDic[mz, charge][1]
            if intens >= intensCutoff:
                if len(finalPepFeatDic) == 0:
                    finalPepFeatDic[mz, charge] = currChgedIsoDic[(mz, charge)]
                else:
                    currFanalList = [x[0] for x in list(finalPepFeatDic.keys()) if x[1] == charge]
                    findBool, matchedMZ = findOneMZ1(mz, currFanalList, 10)
                    if not findBool:
                        finalPepFeatDic[mz, charge] = currChgedIsoDic[(mz, charge)]
                    else:
                        if intens > finalPepFeatDic[matchedMZ, charge][1]:
                            del finalPepFeatDic[matchedMZ, charge]
                            finalPepFeatDic[mz, charge] = currChgedIsoDic[(mz, charge)]

    return finalPepFeatDic


# 将MS1信息写入dic
def trans_ms1to_dic(openedMS1):
    ms1_elute_dic = {}
    i = 0
    f = openedMS1
    while i < len(openedMS1):
        if openedMS1[i][0] != "S":
            i += 1
        else:
            currentScan = int(openedMS1[i].split("\t")[1])
            currentRT = float(openedMS1[i+1].strip().split("\t")[-1])
            currentMzIntens_list = []
            # p = i + 5
            p = i + 1
            while p < len(f):
                if f[p][0] == "S":
                    break
                else:
                    if not f[p][0].isdigit():
                        p += 1
                    else:
                        currentMzIntens_list.append(f[p].rstrip().split(" ")[:2])
                        p += 1
            ms1_elute_dic[currentScan] = [currentRT, currentMzIntens_list]
            i = p
    
    return ms1_elute_dic

# 计算面积
def cal_area(elute_point):
    area = 0
    if len(elute_point) == 1:
        return 0
    else:
        for i in range(len(elute_point) - 1):
            cur_rt, cur_ints = elute_point[i][:2]
            next_rt, next_ints = elute_point[i+1][:2]
            cur_area = (cur_ints + next_ints) * (next_rt - cur_rt) / 2
            area += cur_area
        return round(area, 2)

#寻找对应MS1谱图中的目标离子
def find_mz_in_spec(mz, spec, ttl = 10):
    mz_spec = [float(x[0]) for x in spec]
    deltas = [abs(x- mz) for x in mz_spec]
    min_index = deltas.index(min(deltas))
    if min(deltas) / mz * 1000000 < ttl:
        return True, mz_spec[min_index], float(spec[min_index][1])
    else:
        return False, 0, 0


def elute_curve(mz, charge, RT_start, RT_end, ms1_dic):
    elute_points = []
    for scan in ms1_dic:
        currentRT, currentMzIntens_list = ms1_dic[scan]
        if currentRT >= RT_start  and currentRT <= RT_end :
            sub_spec = [x for x in currentMzIntens_list if float(x[0]) > mz - 4 and  float(x[0]) < mz + 7 ]
            # print(sub_spec)
            if sub_spec == []:
                elute_points.append((currentRT, 0, scan))
            else:
                features, non_iso_spec = detect_isotopic_cluster(sub_spec)
                feated_mzs = [x[0] for x in features if x[1] == charge]
                find_iso_bool, match_mz = findOneMZ1(mz, feated_mzs)
                if find_iso_bool:
                    match_ints = [float(x[1]) for x in sub_spec if float(x[0]) == match_mz][0]
                    elute_points.append((currentRT, match_ints, scan))
                else:
                    if non_iso_spec == []:
                        elute_points.append((currentRT, 0, scan))
                    else:
                        find_bool, match_non_mz, matche_non_intens = find_mz_in_spec(mz, non_iso_spec)
                        if not find_bool:
                            elute_points.append((currentRT, 0, scan))
                        else:
                            elute_points.append((currentRT, matche_non_intens, scan))
    return elute_points


# 主流程
def main():
    # peakbase = 0.1
    tgt = open(tgt_ion_file).readlines()
    for fl in os.listdir(os.getcwd()):
        if fl[-4:] == ".ms1":
            print("the current ms file is " + fl)
            f = open(fl, 'r').readlines()
            repName = fl[:-4] + "_report.csv"
            b = open(repName, 'w')
            b.write(tgt[0][:-1]+",area,scans\n")
            ms1_info_dic = trans_ms1to_dic(f)
            for line in tgt[1:]:
                print(line)
                line_list = line.strip().split(",")
                mz = float(line_list[1])
                charge = int("".join([x for x in line_list[2] if x.isdigit()]))
                RT_start = int("".join([x for x in line_list[3] if x.isdigit()]))
                RT_end = int("".join([x for x in line_list[4] if x.isdigit()]))
                elute_points = elute_curve(mz, charge, RT_start, RT_end, ms1_info_dic)
                area_elute = cal_area(elute_points)
                scans = ";".join([str(x[-1]) for x in elute_points if x[1] != 0])
                b.write("%s,%f,%s\n" % (line.strip(), area_elute, scans))
            b.close()


if __name__ == "__main__":
    main()
    print("Well Done")