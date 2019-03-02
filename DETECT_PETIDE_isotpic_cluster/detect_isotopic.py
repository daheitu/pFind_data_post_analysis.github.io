import os
from copy import deepcopy
from read_spec_from_mgf import read_spec
from delta_ppm import generate_mass_range


mgf_info = read_spec("./test_spec.mgf")
print(mgf_info)
b = open("report.txt", 'w')
cg_delta_dic = {1: 1, 2: 0.5, 3: 0.33333, 4: 0.25, 5: 0.2}
delta_cg_dic = {1:1, 0.5: 2, 0.3333: 3, 0.25:4, 0.2:5}
chr_list = [1, 2, 3, 4, 5]
delta_iso_list = list(cg_delta_dic.values())


def find_cert_mz(mz, ms2_ms_list, begin_idx, tol_ppm):
    k = begin_idx + 1
    [low_ms, up_ms] = generate_mass_range(mz, tol_ppm)
    if ms2_ms_list[-1] < low_ms:
        return False, -1
    else:
        while k < len(ms2_ms_list):
            if float(ms2_ms_list[k]) > up_ms:
                return False, -1
            elif float(ms2_ms_list[k]) >= low_ms:
                return True, k
            else:
                k += 1
# a = [546.5873, 584.9302, 623.2737, 686.6615, 717.3005, 722.9763, 728.9789, 760.6715, 762.0266, 842.0547]
# print(find_cert_mz(762.0266+32, a, 8, 30))


def detect_iso(ms2_dic):
    ms2_info_list = sorted(ms2_dic.items(), key=lambda d:d[1], reverse = True)
    #print(ms2_info_list)
    ms2_ms_list = list(ms2_dic.keys())
    ms2_ms_list.sort()
    chged_iso_dic = {}
    i = 0
    while i < len(ms2_info_list):
        # print(i)
        mz = ms2_info_list[i][0]
        mz_idx = ms2_ms_list.index(mz)
        idx = mz_idx + 1
        b.write("\t".join([str(mz), str(mz_idx)])+"\n")
        if mz_idx > len(ms2_ms_list)-2:
            i += 1
        else:
            find_bool = False
            while idx <= len(ms2_ms_list)-2:
                if float(ms2_ms_list[idx]) > generate_mass_range(float(mz)+1, 20)[1]:
                    break
                else:
                    delta =  float(ms2_ms_list[idx]) - float(mz)
                    delta_min_ther = [abs(x -delta) for x in delta_iso_list]
                    min_delta_idx = delta_min_ther.index(min(delta_min_ther))
                    ther_plus1_mz = float(mz) + delta_iso_list[min_delta_idx]
                    [low_ms, up_ms] = generate_mass_range(ther_plus1_mz, 20)
                    if ms2_ms_list[idx] > low_ms and ms2_ms_list[idx] < up_ms:                    
                        chrg = chr_list[min_delta_idx]
                        interval = delta_iso_list[min_delta_idx]
                        [match_bool, pos] = find_cert_mz(float(mz)+interval*2, ms2_ms_list, idx, 20)
                        if match_bool:
                            find_bool = True
                            iso_cluster = [mz, ms2_ms_list[idx], ms2_ms_list[pos]]
                            m = 3
                            while m < 6:
                                [match_bool, pos] = find_cert_mz(float(mz)+interval*m, ms2_ms_list, idx, 20)
                                if match_bool:
                                    iso_cluster.append(ms2_ms_list[pos])
                                    m += 1
                                else:
                                    break
                            if mz * chrg > 1600:
                                for n in range(1,4):
                                    [match_bool, pos] = find_cert_mz(float(mz)-interval*n, ms2_ms_list, idx-10, 20)
                                    if match_bool:
                                        if ms2_dic[ms2_ms_list[pos]]/ms2_dic[iso_cluster[0]] > 0.3:
                                            iso_cluster.insert(0, ms2_ms_list[pos])
                                        else:
                                            break
                                    else:
                                        break
                            else:
                                pass
                            
                            chged_iso_dic[iso_cluster[0]] = [chrg, iso_cluster]
                            for peak in iso_cluster:
                                ms2_ms_list.remove(peak)
                                ms2_info_list.remove((peak, ms2_dic[peak]))
                            # print(ms2_ms_list)
                            # print(ms2_info_list)
                            b.writelines([str(ele) for ele in ms2_info_list])
                        else:
                            idx += 1
                    else:
                        idx += 1
            if find_bool == True:
                i = i
            else:
                i += 1 
    b.writelines([str(ele) for ele in ms2_info_list])                    
    return chged_iso_dic


def pair_with_Delta_mass(chged_iso_dic, tol_ppm, delta_mass):
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

            [fd_bool, pos] = find_cert_mz(ther_pair_mz, bms_list, i, tol_ppm)
            if fd_bool:
                if chrg not in ft_pair_dic:
                    ft_pair_dic[chrg] = [(bms_list[i], bms_list[pos])]
                else:
                    ft_pair_dic[chrg].append((bms_list[i], bms_list[pos]))
                mz_match_ion = bms_list[pos]
                bms_list.remove(mz)
                bms_list.remove(mz_match_ion)
                i = i
            else:
                i += 1
        if chrg not in chr_ms_left_dic:
            chr_ms_left_dic[chrg] = bms_list
        else:
            print("wrong")
    return ft_pair_dic, chr_ms_left_dic


for ttl in mgf_info:
    ms2_dic = mgf_info[ttl]
    chged_iso_dic = detect_iso(ms2_dic)
    ft_pair_dic, chr_ms_left_dic = pair_with_Delta_mass(chged_iso_dic, 30, 32)
    print(ft_pair_dic)
    print(chr_ms_left_dic)
