

from xml.sax.handler import feature_external_ges


chrgDmassDic = {1: 1, 2: 0.5, 3: 0.33333, 4: 0.25, 5: 0.2}


def is_isotop(mz, real_mz):
    charged_delta_list = list(chrgDmassDic.values())
    delta_mz = abs(mz-real_mz)  #计算质量差绝对值
    delta_list = [abs(x-delta_mz) for x in charged_delta_list]
    min_index = delta_list.index(min(delta_list))
    possible_detla = charged_delta_list[min_index]
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
    intns_mz = [float(x[1]) for x in ms2_spec_list if float(x[0]) == mz][0]
    intns_mz_adj_1 = [float(x[1]) for x in ms2_spec_list if float(x[0]) == mz_adj_1][0]
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
            poss_isoclusters.append([iso_cluter_total, pos_cahrge])
        if len(poss_isoclusters) == 1:
            return poss_isoclusters
        else:
            pos_cls_sorted = sorted(poss_isoclusters, key = lambda x:len(x[0]), reverse= True)
            if len(pos_cls_sorted[1][0]) == len(pos_cls_sorted[0][0]):
                return pos_cls_sorted[:2]
            else:
                return [pos_cls_sorted[0]]
    else:
        bk_pos_charges_mz = look_back(mz, mz_list)
        if bk_pos_charges_mz == []:
            return None
        else:
            poss_isoclusters = []
            for (pos_cahrge, adj_mz) in bk_pos_charges_mz:
                if is_peak_minisOne_reasonable(mz, adj_mz, ms2_spec_list, 0.8):
                    poss_isoclusters.append([[adj_mz, mz], pos_cahrge])
                else:
                    return None



def delet_element(spec_list_ints_order, isocluster_total):
    return [x for x in spec_list_ints_order if float(x[0]) not in isocluster_total]


def detect_isotopic_cluster(ms2_spec_list):
    chr_list = [1, 2, 3, 4, 5]
    mz_list = [float(x[0]) for x in ms2_spec_list]
    spec_list_ints_order = sorted(ms2_spec_list, key=lambda x:float(x[1]), reverse= True) # 按照强度排序
    features = []
    i = 0
    while i < len(spec_list_ints_order):
        mz, ints = [float(x) for x in spec_list_ints_order[i]]
        # print(mz)
        cluster_list = look_and_find_culster(mz, mz_list, ms2_spec_list)
        # print(cluster_list)
        if cluster_list != None:
            for iso_info in cluster_list:
                # print(iso_info)
                iso_cluter_total, pos_cahrge = iso_info[:]
                mono_mz = iso_cluter_total[0]
                most_ints = ints
                features.append((mono_mz, pos_cahrge, mz, most_ints, iso_cluter_total))
                # print(features)
                spec_list_ints_order = delet_element(spec_list_ints_order, iso_cluter_total)
        else:
            i += 1
    return features, spec_list_ints_order


def main():
    f = open('./example_ms1.ms1').readlines()
    i = 0
    while i < len(f):
        if f[i][0] != "S":
            i += 1
        else:
            scan = int(f[i].split("\t")[1])
            ms_spec_list = []
            p = i + 1
            while p < len(f):
                if f[p][0] == "S":
                    break
                else:
                    if not f[p][0].isdigit():
                        p += 1
                    else:
                        ms_spec_list.append(f[p].rstrip().split(" "))
                        p += 1
            features = detect_isotopic_cluster(ms_spec_list)
            print(features)
            i = p


# main()