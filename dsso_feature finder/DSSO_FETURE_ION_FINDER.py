import os

os.chdir(r"G:\DSSO_1008\CID_BASED_Methods\pparse")
relt_ins_thrhd = 0.01
abs_ins_thrhd = 10000
relt_pair_thrhd = 0.3

def generate_ion_mass_range(num, tol):
    deta = num * tol / 1000000
    return [num - deta, num + deta]


def cal_delta_ppm(num1, num2):
    delta = abs(num1- num2)
    avg = (num1 + num2)/2
    ppm = (delta / avg) * 100000
    return ppm 


def cal_pair_rel_int(ins1, ins2):
    ins1 = float(ins1)
    ins2 = float(ins2)
    pair_rel_int = min([ins1, ins2])/max([ins1, ins2])
    return round(pair_rel_int, 2)


def deal_ms2_dic(ms2_dic, delta_mass, ms2_total_dic):
    rep_pair_list = []
    intsy_list = list(ms2_total_dic.values())
    intsy_sorted_list = sorted(intsy_list)
    max_intsy = max(intsy_sorted_list)
    for chrg in ms2_dic:
        delta_mz = delta_mass / int(chrg)
        mz_list = ms2_dic[chrg]

        if len(mz_list) > 1:
            for i in range(len(mz_list) - 1):
                mz = mz_list[i]
                [min_mz, mx_mz] = generate_ion_mass_range(
                    float(mz) + delta_mz, ms2_tol)
                for j in range(i + 1, len(mz_list)):
                    if float(mz_list[j]) > mx_mz:
                        break
                    elif float(mz_list[j]) >= min_mz:
                        mzi_tuple = (float(mz_list[i]),
                                     round(ms2_total_dic[mz_list[i]], 1),
                                     round(
                                             float(ms2_total_dic[mz_list[i]]) /
                                             float(max_intsy), 4),
                                     intsy_sorted_list.index(
                                             ms2_total_dic[mz_list[i]]))
                        mzj_tuple = (float(mz_list[j]),
                                     round(ms2_total_dic[mz_list[j]], 1),
                                     round(float(ms2_total_dic[mz_list[j]]) /
                                            float(max_intsy), 4),
                                     intsy_sorted_list.index(
                                             ms2_total_dic[mz_list[j]]))
                        rep_pair_list.append(
                            '[' + str(mzi_tuple) + ', ' + str(mzj_tuple) + ']')
                    else:
                        pass
        else:
            continue
    return rep_pair_list


# f = open("BSA_DSSO_BIO_R1_T1_HCDFT_deisotopic.mgf", 'r').readlines()
# b = open("DSSO_rep.txt", "w")
ms2_tol = 30

file_list = os.listdir(os.getcwd())
for fl in file_list:
    if fl[-14:] == "deisotopic.mgf":
        name = fl[:-14] + "rep_pair.txt"
        f = open(fl, 'r').readlines()
        b = open(name, 'w')

        i = 0
        while i < len(f):
            ms2_total_dic = {}
            ms2_dic = {}
            if f[i].strip() == "BEGIN IONS":
                title = f[i + 1].strip().split("=")[1]
            else:
                print("wrong")
            p = i + 4
            if title[-5:] == "0.dta":
                precusor_charge = int(f[i + 2].strip().split("=")[1][:-1])
                pepMZ = float(f[i + 3].strip().split("=")[1])
                pepMass = round((pepMZ - 1.00782) * precusor_charge, 4)

                while p < len(f) and f[p][0].isdigit():
                    line_list = f[p].strip().split(" ")
                    p += 1
                    if len(line_list) == 2:
                        ms2_total_dic[line_list[0]] = float(line_list[1])
                    else:
                        chrg = line_list[-1]
                        mz = line_list[0]
                        its = line_list[1]
                        ms2_total_dic[mz] = float(its)
                        if chrg not in ms2_dic:
                            ms2_dic[chrg] = [mz]
                        else:
                            ms2_dic[chrg].append(mz)
                rep_pair_list = deal_ms2_dic(ms2_dic, 32, ms2_total_dic)
                if rep_pair_list:
                    wt_list = [
                        title,
                        str(precusor_charge),
                        str(pepMZ),
                        str(pepMass), 
                        ";".join(rep_pair_list)
                    ]
                    # print(rep_pair_list)
                else:
                    wt_list = [
                        title,
                        str(precusor_charge),
                        str(pepMZ),
                        str(pepMass)
                    ]

                b.write("\t".join(wt_list))
                b.write("\n")
            else:
                while p < len(f) and f[p][0].isdigit():
                    p += 1
            i = p + 1

        b.close()
        pair_num_dic = {}
        pair_num_rlInt_dic = {}
        pair_num_rlPint_dic = {}
        pair_num_ppm_dic = {}
        pair_num_IntRank_dic = {}
        b = open(name, 'r').readlines()
        print(b[1])
        pair_num_dic[0] = 0
        for line in b:
            line_list = line.strip().split("\t")
            if len(line_list) == 4:
                pair_num_dic[0] += 1
            else:
                pair_list = line_list[-1].split(";")
                pair_list_len = len(pair_list)
                if pair_list_len not in pair_num_dic:
                    pair_num_dic[pair_list_len] = 1
                else:
                    pair_num_dic[pair_list_len] += 1
                for pair in pair_list:
                    mz1_info_list = pair[1:-1].split("), (")[0][1:].split(",")
                    mz2_info_list = pair[1:-1].split("), (")[1][:-1].split(",")
                    delta_ppm = cal_delta_ppm(float(mz1_info_list[0]), float(mz2_info_list[0]))
                    pair_rl_int = cal_pair_rel_int(mz1_info_list[1], mz2_info_list[1])
                    rlInt_list = [float(mz1_info_list[2]), float(mz2_info_list[2])]
                    rank_delta = abs(int(mz1_info_list[-1]) - int(mz2_info_list[-1]))
                    if pair_list_len not in pair_num_ppm_dic:
                        pair_num_ppm_dic[pair_list_len] = [delta_ppm]
                    else:
                        pair_num_ppm_dic[pair_list_len].append(delta_ppm)
                    
                    if pair_list_len not in pair_num_rlPint_dic:
                        pair_num_rlPint_dic[pair_list_len] = [pair_rl_int]
                    else:
                        pair_num_rlPint_dic[pair_list_len].append(pair_rl_int)
                    
                    if pair_list_len not in pair_num_rlInt_dic:
                        pair_num_rlInt_dic[pair_list_len] = rlInt_list
                    else:
                        pair_num_rlInt_dic[pair_list_len].extend(rlInt_list)

                    if pair_list_len not in pair_num_IntRank_dic:
                        pair_num_IntRank_dic[pair_list_len] = [rank_delta]
                    else:
                        pair_num_IntRank_dic[pair_list_len].append(rank_delta)

        print(pair_num_dic)
        print(pair_num_rlInt_dic)
        print(pair_num_rlPint_dic)
        print(pair_num_ppm_dic)
        print(pair_num_IntRank_dic)
        """
        stac = open(fl[:-14] + "stac.txt", 'w')
        pair_num_list = sorted(list(pair_num_dic.keys()))
        for pair_num in pair_num_list:
            stac.write("\t".join([str(pair_num),
                                  str(pair_num_dic[pair_num])]) + "\n")

        stac.close()
        """