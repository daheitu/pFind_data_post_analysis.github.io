from ast import Return
from operator import truediv
import os
from copy import deepcopy

from urllib3 import Retry
from read_spec_from_mgf import read_spec
from delta_ppm import generate_mass_range

#os.chdir(r"D:\Onedriver\OneDrive\github\pFind_data_post_analysis.github.io\DETECT_PETIDE_isotpic_cluster")

chrgDmassDic = {1: 1, 2: 0.5, 3: 0.33333, 4: 0.25, 5: 0.2}

def compareTwoNum(realNum, TheroNum, tolPPM):
    if abs(TheroNum- realNum)/TheroNum * 1000000 < tolPPM:
        return True
    else:
        return False


def findOneMZ(mz, ms2MzList, begin_idx, tol_ppm):
    low_ms, up_ms = generate_mass_range(mz, tol_ppm)
    #
    if mz >  ms2MzList[begin_idx]:
        k = begin_idx + 1
        
        if ms2MzList[-1] < low_ms:
            return False, -1
        else:
            while k < len(ms2MzList):
                if float(ms2MzList[k]) > up_ms:
                    return False, -1
                elif float(ms2MzList[k]) >= low_ms:
                    return True, ms2MzList[k]
                else:
                    k += 1
            if k == len(ms2MzList):
                return False, -1
    else:
        k = begin_idx - 1
        if ms2MzList[0] > up_ms:
            return False, -1
        else:
            while k > 0:
                if float(ms2MzList[k]) < low_ms:
                    return False, -1
                elif float(ms2MzList[k]) <= up_ms:
                    return True, ms2MzList[k]
                else:
                    k -= 1
            if k == -1:
                return False, -1


def is_isotop(mz, real_mz):
    charged_delta_list = list(chrgDmassDic.values())
    delta_mz = abs(mz-real_mz)  #计算质量差绝对值
    delta_list = [abs(x-delta_mz) for x in charged_delta_list]
    min_index = delta_list.index(min(delta_list))
    possible_detla = charged_delta_list(min_index)
    possible_charge = int(1.00/ possible_detla)
    if mz > real_mz:
        thero_mz = mz - possible_detla
        if abs((thero_mz - real_mz) / thero_mz) < 1e-5:
            return True, possible_charge
        else:
            return False, -1
    else:
        thero_mz = mz + possible_detla
        if abs((thero_mz - real_mz) / thero_mz) < 1e-5:
            return True, possible_charge
        else:
            return False, -1

# 在特定窗口里寻找mass
def find_tgt_mass(tgt_mz, mz_list, mass_start, mass_window, direction):
    if direction == "back":
        mass_window = [mass_start - mass_window, mass_start]
    elif direction == "forwad":
        mass_window = [mass_start, mass_start + mass_window]
    
    for mass in mz_list:
        if mass > mass_window[0] and mass < mass_window[1]:
            ppm_delta = (mass - tgt_mz) / tgt_mz
            if abs(ppm_delta) < 1e-5:
                return True, mass
            else:
                continue
    return False, -1
    




# 从小窗口寻找其同位素峰
def look_for_clster(mz, adj_mz, charge, mz_list):
    if mz > adj_mz:
        iso_cluster = [adj_mz, mz]
        mass_window = 4
        delta = 1 / charge
        for i in range(1, 3):
            cur_iso_mz = adj_mz - delta * i
            find_bool, iso_found = find_tgt_mass(cur_iso_mz, mz_list, mz, mass_window, "back")
            if find_bool:
                iso_cluster.insert(0, iso_found)
            else:
                break
        return iso_cluster
    else:
        iso_cluster = [mz, adj_mz]
        mass_window = 7
        delta = 1 / charge
        for i in range(1, 6):
            cur_iso_mz = adj_mz + delta * i
            find_bool, iso_found = find_tgt_mass(cur_iso_mz, mz_list, mz, mass_window, "forwad")
            if find_bool:
                iso_cluster.append(iso_found)
            else:
                break
        return iso_cluster



# 正向寻找同位素峰
def look_forwad(mz, mz_list):
    pos_charges_mz = []
    for mass in mz_list:
        if mass > mz and mass < mz + 1.1:
            find_bool, poscharge = is_isotop(mz, mass)
            if find_bool:
                pos_charges_mz.append((poscharge, mass))
    return pos_charges_mz


def look_back(mz, mz_list):
    pos_charges_mz = []
    for mass in mz_list:
        if mass < mz and mass > mz - 1.1:
            find_bool, poscharge = is_isotop(mz, mass)
            if find_bool and poscharge == 1:
                pos_charges_mz.append((poscharge, mass))
    return pos_charges_mz


def is_peak_minisOne_reasonable(mz, mz_adj_1, ms2_spec_list, cut_off=0.4):
    intns_mz = [float(x[1]) for x in ms2_spec_list if float(x) == mz]
    intns_mz_adj_1 = [float(x[1]) for x in ms2_spec_list if float(x) == mz_adj_1]
    if intns_mz_adj_1 / intns_mz > cut_off:
        return True
    else:
        return False



def look_and_find_culster(mz, mz_list, ms2_spec_list, tol = 10):
    fw_pos_charges_mz = look_forwad(mz, mz_list)
    if fw_pos_charges_mz != []:
        poss_isoclusters = []
        for (pos_cahrge, adj_mz) in fw_pos_charges_mz:
            iso_cluster_fw = look_for_clster(mz, adj_mz, pos_cahrge, mz_list)
            iso_cluster_bk = look_for_clster(adj_mz, mz, pos_cahrge, mz_list)
            if len(iso_cluster_bk) > 2:
                if is_peak_minisOne_reasonable(mz, iso_cluster_bk[-3], ms2_spec_list):
                    iso_cluter_total = iso_cluster_bk[:-2]
                    iso_cluter_total.extend(iso_cluster_fw)
                else:
                    iso_cluter_total = iso_cluster_fw
            else:
                iso_cluter_total = iso_cluster_fw
            # mono_mz = iso_cluter_total[0]
            # num_iso_cluster_peaks = len(iso_cluter_total)
            poss_isoclusters.append((iso_cluter_total, pos_cahrge))
        if len(poss_isoclusters) == 1:
            return poss_isoclusters
        else:
            pos_cls_sorted = sorted(poss_isoclusters, key = lambda x:len(x[0]), reverse= True)
            if len(pos_cls_sorted[1][0]) == len(pos_cls_sorted[0][0]):
                return pos_cls_sorted[:2]
            else:
                return pos_cls_sorted[0]
    else:
        bk_pos_charges_mz = look_back(mz, mz_list)
        if bk_pos_charges_mz == []:
            return None
        else:
            poss_isoclusters = []
            for (pos_cahrge, adj_mz) in bk_pos_charges_mz:
                if is_peak_minisOne_reasonable(mz, adj_mz, ms2_spec_list, 0.8):
                    poss_isoclusters.append(([adj_mz, mz], pos_cahrge))
                else:
                    return None
            





def detect_isotopic_cluster(ms2_spec_list):
    chr_list = [1, 2, 3, 4, 5]
    mz_list = [float(x[0]) for x in ms2_spec_list]
    spec_list_ints_order = sorted(ms2_spec_list, key=lambda x:float(x[1]), reverse= True) # 按照强度排序
    features = []
    i = 0
    while i < len(spec_list_ints_order):
        mz, ints = [float(x) for x in spec_list_ints_order[i]]
        cluster_list = look_and_find_culster(mz, mz_list, ms2_spec_list)
        if cluster_list != None:
            for (iso_cluter_total, pos_cahrge) in cluster_list:
                mono_mz = iso_cluter_total[0]
                features.append((mono_mz, pos_cahrge, mz, iso_cluter_total))
                for mz_c in iso_cluter_total:
                    if mz_c in [x[0] for x in spec_list_ints_order]:
                        tgt_ele = [x for x in spec_list_ints_order if float(x[0]) == mz_c]
                        spec_list_ints_order.remove(tgt_ele)
        else:
            i += 1
    return features




# 去同位素峰
def detectIsotopic(ms2_dic):
    chrgDmassDic = {1: 1, 2: 0.5, 3: 0.33333, 4: 0.25, 5: 0.2}
    chr_list = [1, 2, 3, 4, 5]
    dmassChargeList = list(chrgDmassDic.values())
    ms2_info_list = sorted(ms2_dic.items(), key=lambda d:d[1], reverse = True)
    ms2MzList = list(ms2_dic.keys())
    ms2MzList.sort()
    chged_iso_dic = {}
    i = 0
    while i < len(ms2_info_list):
        mz = ms2_info_list[i][0]
        mz_idx = ms2MzList.index(mz)
        idx = mz_idx + 1
        # fileToWrite.write("\t".join([str(mz), str(mz_idx)])+"\n")
        if mz_idx > len(ms2MzList)-2:
            i += 1
        else:
            find_bool = False
            while idx <= len(ms2MzList)-2:
                if float(ms2MzList[idx]) > generate_mass_range(float(mz)+1, 20)[1]:
                    break
                else:
                    delta =  float(ms2MzList[idx]) - float(mz)
                    delta_min_ther = [abs(x - delta) for x in dmassChargeList]
                    min_delta_idx = delta_min_ther.index(min(delta_min_ther))
                    ther_plus1_mz = float(mz) + dmassChargeList[min_delta_idx]
                    mzPlus1Bool = compareTwoNum(ms2MzList[idx], ther_plus1_mz, 20)
                    intensPlus1Bool = ms2_dic[ms2MzList[idx]]/ms2_dic[mz] > 0.3
                    if mzPlus1Bool and intensPlus1Bool:
                        chrg = chr_list[min_delta_idx]
                        interval = dmassChargeList[min_delta_idx]
                        #print(interval)
                        match_bool, matchedMZ = findOneMZ(float(mz)+interval*2, ms2MzList, idx, 20)
                        if match_bool:
                            find_bool = True
                            iso_cluster = [mz, ms2MzList[idx], matchedMZ]
                            # print(iso_cluster)
                            m = 3
                            while m < 6:
                                [matchMoreBool, matchMoreMZ] = findOneMZ(float(mz)+interval*m, ms2MzList, idx, 20)
                                if matchMoreBool:
                                    iso_cluster.append(matchMoreMZ)
                                    m += 1
                                else:
                                    break
                            if iso_cluster[0] * chrg > 1600:                                
                                for n in range(1,4):
                                    [matchLessBool, matchLessMZ] = findOneMZ(iso_cluster[0]-interval, ms2MzList, idx, 20)
                                    #print(matchLessBool, matchLessMZ)
                                    if matchLessBool:
                                        if ms2_dic[matchLessMZ]/ms2_dic[iso_cluster[0]] > 0.3:
                                            iso_cluster.insert(0, matchLessMZ)
                                        else:
                                            break
                                    else:
                                        break
                            else:
                                pass
                            # print(iso_cluster)
                            chged_iso_dic[iso_cluster[0]] = [chrg, iso_cluster]
                            for peak in iso_cluster:
                                if peak in ms2MzList:
                                    ms2MzList.remove(peak)
                                    ms2_info_list.remove((peak, ms2_dic[peak]))
                                else:
                                    pass
                            # fileToWrite.writelines([str(ele) for ele in ms2_info_list])
                        else:
                            idx += 1
                    else:
                        idx += 1
            if find_bool == True:
                i = i
            else:
                i += 1 
    # fileToWrite.writelines([str(ele) for ele in ms2_info_list])                    
    return chged_iso_dic


def pairedwithDeltaMass(chged_iso_dic, tol_ppm, delta_mass):
    ft_pair_dic = {}
    chr_ms_left_dic = {}
    chr_ms_dic = {}
    for key in chged_iso_dic:
        charge = chged_iso_dic[key][0]
        if charge not in chr_ms_dic:
            chr_ms_dic[charge] = [key]
        else:
            chr_ms_dic[charge].append(key)
    for chrg in chr_ms_dic:
        chr_ms_dic[chrg].sort()
        delta_mz = delta_mass/chrg
        bms_list = deepcopy(chr_ms_dic[chrg])
        
        i = 0        
        while i < len(bms_list)-1:
            mz = bms_list[i]
            ther_pair_mz = mz + delta_mz

            [fd_bool, realMZ] = findOneMZ(ther_pair_mz, bms_list, i, tol_ppm)
            if fd_bool:
                if chrg == 1:
                    [fd50_bool, real50MZ] = findOneMZ(ther_pair_mz+18.010, bms_list, i+1, tol_ppm)
                    if fd50_bool:
                        if chrg not in ft_pair_dic:
                            ft_pair_dic[chrg] = [(bms_list[i], realMZ, real50MZ)]
                        else:
                            ft_pair_dic[chrg].append((bms_list[i],realMZ, real50MZ))
                        bms_list.remove(mz)
                        bms_list.remove(realMZ)
                        bms_list.remove(real50MZ)
                    else:
                        if chrg not in ft_pair_dic:
                            ft_pair_dic[chrg] = [(bms_list[i],realMZ)]
                        else:
                            ft_pair_dic[chrg].append((bms_list[i],realMZ))
                        bms_list.remove(mz)
                        bms_list.remove(realMZ)
                else:        
                    if chrg not in ft_pair_dic:
                        ft_pair_dic[chrg] = [(bms_list[i], realMZ)]
                    else:
                        ft_pair_dic[chrg].append((bms_list[i],realMZ))
                    bms_list.remove(mz)
                    bms_list.remove(realMZ)
                i = i
            else:
                i += 1
        if chrg not in chr_ms_left_dic:
            chr_ms_left_dic[chrg] = bms_list
        else:
            print("wrong")
    return ft_pair_dic, chr_ms_left_dic


def statistcPairNum(ft_pair_dic):
    print(ft_pair_dic)
    pairNum = 0
    for chrg in ft_pair_dic:
        pairNum += len(ft_pair_dic[chrg])
    
    return pairNum


def deal_deisotopiced_mgf(mgfFl):
    ft_pair_dic ={}
    if mgfFl[-14:] != "deisotopic.mgf":
        print('wrong')
        return ft_pair_dic
    else:
        fl = open(mgfFl, 'r').readlines()
        
        i = 0
        while i < len(fl):
            if fl[i].strip() != "BEGIN IONS":
                i += 1
            else:
                scan = fl[i+1].split(".")[1]
                pairDic = {}
                p = i + 5
                while p < len(fl):
                    if fl[p].strip() == "END IONS":
                        break
                    else:
                        lineList = fl[p].strip().split(' ')
                        if len(lineList) == 3:
                            chrg = int(lineList[-1])
                            if chrg not in pairDic:
                                pairDic[chrg] = [float(lineList[0])]
                            else:
                                pairDic[chrg].append(float(lineList[0]))
                        p += 1
                print(pairDic)
                pickedPairDic = pairedwithDeltaMass(pairDic, 30, 32)[0]
                pairedNum  = statistcPairNum(pickedPairDic)
                #print(pairedNum)
                ft_pair_dic[scan] = pairedNum
                i = p+1
    return ft_pair_dic


# print(deal_deisotopiced_mgf('rep_dsso_deisotopic.mgf'))


def main():
    mgf_info = read_spec("./test1.mgf")
    # print(mgf_info)
    b = open("report.txt", 'w')
    pairNumDic = {}
    for ttl in mgf_info:
        ms2_dic = mgf_info[ttl]
        chged_iso_dic = detectIsotopic(ms2_dic)
        print(chged_iso_dic)
        t_pair_dic, chr_ms_left_dic = pairedwithDeltaMass(chged_iso_dic, 30, 31.97)
        # print(t_pair_dic)
        pairNum = statistcPairNum(t_pair_dic)
        if pairNum not in pairNumDic:
            pairNumDic[pairNum] = 1
        else:
            pairNumDic[pairNum] += 1
    for num in sorted(list(pairNumDic.keys())):
        b.write(str(num) + '\t' + str(pairNumDic[num]) + '\n')
    b.close()


if __name__ == "__main__":
    main()
