import os

os.chdir(r"G:\DSSO_1008\CID_BASED_Methods\pparse")

f = open("DSSO_CID_MS2_CID_MS3_T1_CIDFT_deisotopic.mgf", 'r').readlines()
b = open("DSSO_rep.txt", "w")
ms2_tol = 30


def generate_ion_mass_range(num, tol):
    deta = num * tol / 1000000
    return [num - deta, num + deta]


def deal_ms2_dic(ms2_dic, delta_mass):
    rep_pair_list = []
    for chrg in ms2_dic:
        delta_mz = delta_mass / int(chrg)
        mz_list = []
        for word in ms2_dic[chrg]:
            mz_list.append(word[0])
        if len(mz_list) > 1:
            for i in range(len(mz_list) - 1):
                mz = mz_list[i]
                [min_mz, mx_mz] = generate_ion_mass_range(
                    float(mz) + delta_mz, ms2_tol)
                for j in range(i + 1, len(mz_list)):
                    if float(mz_list[j]) > mx_mz:
                        break
                    elif float(mz_list[j]) >= min_mz:
                        rep_pair_list.append(str((mz_list[i], mz_list[j])))
        else:
            continue
    return rep_pair_list


i = 0
while i < len(f):
    ms2_dic = {}
    if f[i].strip() == "BEGIN IONS":
        title = f[i + 1].strip().split("=")[1]
        precusor_charge = int(f[i + 2].strip().split("=")[1][:-1])
        pepMZ = float(f[i + 3].strip().split("=")[1])
        pepMass = round((pepMZ - 1.00782) * precusor_charge, 4)
    else:
        print("wrong")
    print(i, title)
    p = i + 4
    while p < len(f) and f[p][0].isdigit():
        line_list = f[p].strip().split(" ")
        p += 1
        if len(line_list) == 2:
            continue
        else:
            chrg = line_list[-1]
            mz = line_list[0]
            its = line_list[1]
            if chrg not in ms2_dic:
                ms2_dic[chrg] = [(mz, its)]
            else:
                ms2_dic[chrg].append((mz, its))
    rep_pair_list = deal_ms2_dic(ms2_dic, 32)
    if rep_pair_list:
        print(rep_pair_list)
    else:
        pass
    wt_list = [title, str(precusor_charge), str(pepMZ), \
        str(pepMass), ",".join(rep_pair_list)]
    b.write("\t".join(wt_list))
    b.write("\n")
    i = p + 1

b.close()