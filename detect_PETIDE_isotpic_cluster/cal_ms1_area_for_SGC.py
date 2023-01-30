# coding:utf-8


import os
import random
from detect_isotopics_20220405 import detect_isotopic_cluster
import gc


gc.enable()

# os.chdir(r'G:\test_mthe')
ms1_path = r"G:\test_sgc\20220303_SP3_SP36_desalt_IP_ProA_DMP_2517_Elution_reject1_4uL.ms1"
rep_xl_pep_path = r"G:\test_sgc\20220303_SP3_SP36_desalt_IP_ProA_DMP_2517_Elution_reject1_4uL_test\reports\Total_Pmt3_con_2022.04.06.filtered_cross-linked_peptides.csv"

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


# 将MS1信息写入dic
def trans_ms1to_dic(openedMS1):
    ms1_elute_dic = {}
    i = 0
    f = openedMS1
    while i < len(openedMS1):
        if openedMS1[i][0] != "S":
            i += 1
        else:
            print(i)
            currentScan = int(f[i].split("\t")[1])
            currentRT = float(f[i+1].strip().split("\t")[-1])
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


# 有改进空间，可以以保留时间为准
def find_mz_in_ms1(mz, charge, scan_star, scan_end, ms1_dic):
    elute_points = []
    scans_ms1 = list(ms1_dic.keys())
    delta_list = [abs(x - scan_star) for x in scans_ms1]
    # min_idx = delta_list.index(min(delta_list))
    for scan in ms1_dic.keys():
        if scan > scan_star and scan < scan_end:
            app_RT, appro_mzintns_list = ms1_dic[scan]
            appro_features = detect_isotopic_cluster(appro_mzintns_list)
            tgt_mz_list = [x[0] for x in appro_features if x[1] == charge]
            find_bool, matched_mz = findOneMZ1(mz, tgt_mz_list, 10)
            if find_bool:
                intns_matched = [x[3] for x in appro_features if x[0] == matched_mz and x[1] == charge][0]
                elute_points.append((app_RT, intns_matched, scan))
            else:
                elute_points.append((app_RT, 0, scan))
    return elute_points


def cal_elute_area(points):
    if len(points) < 2:
        return 0
    else:
        area = 0
        max_intns = max([x[1] for x in points])
        for i in range(len(points) - 1):
            cur_rt, cur_ints = points[i][:2]
            next_rt, next_ints = points[i+1][:2]
            if max([cur_ints, next_ints]) >= max_intns * 0.05:
                cur_area = (cur_ints + next_ints) * (next_rt - cur_rt) / 2
                area += cur_area
        return round(area, 2)


# 从肽段文件读取mz charge信息
def read_xl_pep_file(xl_path, ms1_dic, rep_name):
    f = open(xl_path).readlines()
    b = open(rep_name, 'w')
    b.write(f[0])
    b.write(",,mz,charge,area,scans_MS1, scans_ms2\n")
    
    i = 2
    while i < len(f):
        if not f[i][0].isdigit():
            i += 1
        else:
            b.write(f[i])
            print("the line %d in f" % i)
            pep_mass = float(f[i].split(",")[2])
            ion_dic = {}
            p = i + 1
            while p < len(f):
                if f[p][0].isdigit():
                    break
                else:
                    linelist = f[p].split(",")
                    title_spec = linelist[2]
                    # print(title_spec)
                    mix_idx = title_spec.split(".")[-2]
                    scan = int(title_spec.split(".")[-4])
                    if mix_idx == "0":
                        charge = int(linelist[3])
                        pre_mass = float(linelist[4])
                        pre_mz = round((pre_mass + (charge - 1) * 1.0078) / charge, 5)
                        if ion_dic == {}:
                            ion_dic[pre_mz, charge] = [scan]
                        else:
                            currFanalList = [x[0] for x in list(ion_dic.keys()) if x[1] == charge]
                            findBool, matchedMZ = findOneMZ1(pre_mz, currFanalList, 10)
                            if not findBool:
                                ion_dic[pre_mz, charge] = [scan]
                            else:
                                aver_mz = round((pre_mz + matchedMZ) / 2, 5)
                                scans = ion_dic[matchedMZ, charge]
                                scans.append(scan)
                                del ion_dic[matchedMZ, charge]
                                ion_dic[aver_mz, charge] = scans
                
                    p += 1
            # print(ion_dic)
            for mz,charge in ion_dic:
                # print("the current ion is %f %d" % (mz, charge))
                scans_ms2 = sorted(ion_dic[mz, charge])
                scan_min = min(scans_ms2) - 200
                scan_max = max(scans_ms2) + 200
                elute_points = find_mz_in_ms1(mz, charge, scan_min, scan_max, ms1_dic)
                elute_area = cal_elute_area(elute_points)
                ms1_scans = [x[2] for x in elute_points if x[1] != 0]
                b.write(",,%f,%d,%f,%s,%s\n"% (mz,charge, elute_area,";".join([str(x) for x in ms1_scans]),";".join([str(x) for x in scans_ms2])))

            i = p
    b.close()
    

def main():
    rep_path = os.path.join(os.path.dirname(rep_xl_pep_path), "summary_areas.csv")
    f = open(ms1_path).readlines()
    ms1_dic = trans_ms1to_dic(f)
    read_xl_pep_file(rep_xl_pep_path, ms1_dic, rep_path)
           

if __name__ == "__main__":
    main()
    print("Well Done")