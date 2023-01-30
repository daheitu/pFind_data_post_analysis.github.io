import imp
import os
import random
from detect_isotopics_20220405 import detect_isotopic_cluster
import gc


gc.enable()

os.chdir(r'G:\test_mthe')

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


def detectGlobalBaseLine(ms1_elute_dic, mz, randNum):
    scans_list = sorted(list(ms1_elute_dic.keys()))
    totalScan = scans_list[-1]
    intensList = []
    for i in range(0, randNum):
        seed = random.randint(1, totalScan)
        for scan in scans_list:
            if scan > seed - 200 and scan < seed + 200:
                cur_mz_ints_list = ms1_elute_dic[scan][1]
                cur_mz_list = [float(x[0]) for x in cur_mz_ints_list]
                find_bool, match_mz = findOneMZ1(mz, cur_mz_list)
                if find_bool:
                    intes_matched = [x[1] for x in cur_mz_ints_list if float(x[0]) == match_mz][0]
                    intensList.append(float(intes_matched))
    
    return round(sum(intensList)/len(intensList), 1)


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


# 计算每个肽段洗脱峰面积
def cal_eluteArea(finalPepFeatDic, ms1_dic, scanRange, repName):
    b = open(repName, 'w')
    b.write("\t".join(["mz, charge,  Max_Intensity, max@scanNum, elute_Area@max, scans"]) + "\n")
    scans_list = list(ms1_dic.keys())
    scans_list.sort()
    for mz,charge in finalPepFeatDic:
        print("the current peptide is %f, %d" % (mz, charge))
        max_intn = finalPepFeatDic[mz, charge][1]
        scan_cur = finalPepFeatDic[mz, charge][-1]
        elute_point = []
        fw_stop_point = 0
        for scan in scans_list:
            if fw_stop_point == 2:
                break
            else:
                if scan >= scan_cur and scan < scan_cur + scanRange:
                    RT_cur, mzIntns_list = ms1_dic[scan]
                    currChgedIsoFeatrues = detect_isotopic_cluster(mzIntns_list)
                    fea_mz_list_cur = [x[0] for x in currChgedIsoFeatrues if x[1] == charge]
                    find_bool, matched_mz = findOneMZ1(mz, fea_mz_list_cur)
                    if find_bool:
                        matched_ints = [x[3] for x in currChgedIsoFeatrues if x[1] == charge and x[0] == matched_mz][0]
                        if matched_ints / max_intn > 0.05:
                            elute_point.append((RT_cur, matched_ints, scan))
                        else:
                            elute_point.append((RT_cur, matched_ints, scan))
                            fw_stop_point += 1
                    else:
                        fw_stop_point += 1
        bk_stop_point = 0
        for scan in scans_list[::-1]:
            if bk_stop_point == 2:
                break
            else:
                if scan < scan_cur and scan > scan_cur - scanRange:
                    RT_cur, mzIntns_list = ms1_dic[scan]
                    currChgedIsoFeatrues = detect_isotopic_cluster(mzIntns_list)
                    fea_mz_list_cur = [x[0] for x in currChgedIsoFeatrues if x[1] == charge]
                    find_bool, matched_mz = findOneMZ1(mz, fea_mz_list_cur)
                    if find_bool:
                        matched_ints = [x[3] for x in currChgedIsoFeatrues if x[1] == charge and x[0] == matched_mz][0]
                        if matched_ints / max_intn > 0.05:
                            elute_point.insert(0, (RT_cur, matched_ints, scan))
                        else:
                            elute_point.insert(0, (RT_cur, matched_ints, scan))
                            bk_stop_point += 1
                    else:
                        bk_stop_point += 1
        # print(elute_point)
        area = cal_area(elute_point)
        print(area)
        wlist = [mz, charge, max_intn, scan_cur, area, ";".join([str(x[2]) for x in elute_point])]
        b.write(",".join([str(x) for x in wlist])+ "\n")
    
    b.close()


def main():
    # peakbase = 0.1
    for fl in os.listdir(os.getcwd()):
        if fl[-4:] == ".ms1":
            print("the current ms file is " + fl)
            f = open(fl, 'r').readlines()
            repName = fl[:-4] + "pep_list_perc5.csv"
            ms1_info_dic = trans_ms1to_dic(f)
            baseIntens = detectGlobalBaseLine(ms1_info_dic, 462.14658, 5)
            finalPepFeatDic = detectPepFeat(ms1_info_dic, baseIntens)
            cal_eluteArea(finalPepFeatDic, ms1_info_dic, 200, repName)


if __name__ == "__main__":
    main()