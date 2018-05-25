import os

os.chdir(r"E:\MASS SEPTRA DATA\20170310\H1")
deta_ppm = 5
feature_ion_mass_list = [147.1128, 244.1656, 391.2340]

LK_7_list = [147.1128, 244.1656, 391.2340]
gr_N_list = [175.1189, 246.1561, 317.1932, 388.2303]
gr_C_list = [185.1285, 256.1656, 327.2027, 398.2398]
CR_C_list = [175.1189, 274.1874, 421.2558]
DK_list = [147.1128, 261.1557, 374.2398, 329.1568]


def generate_ion_mass_range(num):
    deta = num * deta_ppm / 100000
    return num - deta, num + deta


def words_frequcecy(list):
    mass_frequcy_dic = {}
    for word in list:
        if word not in mass_frequcy_dic:
            mass_frequcy_dic[word] = 0
        mass_frequcy_dic[word] += 1
    return sorted(mass_frequcy_dic.items(), key=lambda d: d[1], reverse=True)


def judge_spectra(mass_list):
    n = 0
    for ion_mass in feature_ion_mass_list:
        down_mass, up_mass = generate_ion_mass_range(ion_mass)
        bool_feature_ion = False
        for fragm_mass in mass_list:
            if float(fragm_mass) > down_mass and float(fragm_mass) < up_mass:
                bool_feature_ion = True
                break
            else:
                continue

        if bool_feature_ion is False:
            break
        else:
            n += 1
            continue

    if n == len(feature_ion_mass_list):
        return True
    else:
        return False


file_list = os.listdir(os.getcwd())
for name in file_list:
    if name[-3:] == "mgf":
        f = open(name)
        all = f.readlines()
        spec_be_list = []
        spec_end_list = []
        for i in range(len(all)):
            if all[i].strip() == "BEGIN IONS":
                spec_be_list.append(i)
            elif all[i].strip() == "END IONS":
                spec_end_list.append(i)
        print(len(spec_be_list))
        if len(spec_be_list) != len(spec_end_list):
            print("erro")
        else:
            report_file_name = name[:name.find(".mgf")] + ".txt"
            b = open(report_file_name, 'w')
            b.write("\t".join(["title", "charge", "m/z", "mass"]))
            b.write("\n")
            for k in spec_be_list:
                title = all[k + 1].strip()
                charge = all[k + 2][7:9]
                mass_over_z = all[k + 4][all[k + 4].find("=") + 1:all[k + 4]
                                         .find("\n")]
                mass = float(mass_over_z) * int(charge[0]) - int(charge[0])
                w_list = [title, charge, mass_over_z, str(mass)]
                ms2_ion_list = []
                for m in range(k + 5, spec_end_list[spec_be_list.index(k)]):
                    ms2_ion_list.append(all[m].strip().split(" ")[0])

                if judge_spectra(ms2_ion_list):
                    b.write("\t".join(w_list))
                    b.write("\n")
            b.close()

        b = open(report_file_name, 'a+')
        b.seek(0)
        # print(len(b.readlines()))
        mass_list = []
        table = open(report_file_name).readlines()
        for i in range(1, len(table)):
            mass_list.append(round(float(table[i].split("\t")[3]), 1))

        mass_frequcy_dic = words_frequcecy(mass_list)
        print(mass_frequcy_dic)
        for key in mass_frequcy_dic:
            b.write("\t".join([str(key[0]), str(key[1])]))
            b.write("\n")
        b.close()
        f.close()
