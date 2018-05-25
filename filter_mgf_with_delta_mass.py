import os

os.chdir(r"E:\MASS SEPTRA DATA\20170310\H1")
modi_mass_range = [-10, 400]
deta_ppm = 10
noise_cutoff = 0.10
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


def get_intensity_threshhold(num1, num2, all):
    ms2_ion_inten_list = []
    for m in range(num1 + 5, num2):
        ms2_ion_inten_list.append(float(all[m].strip().split(" ")[0]))
    return max(ms2_ion_inten_list) * noise_cutoff


def judge_spectra(feature_list, mass_list):
    n = 0
    for ion_mass in feature_list[1:]:
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

    if n == len(feature_ion_mass_list[1:]):
        return True
    else:
        return False


def get_spec_info(num1, num2, all):
    title = all[num1 + 1].strip()
    charge_state = all[num1 + 2][7:9]
    charge = int(charge_state[0])
    mass_over_z = all[num1 + 4].strip().split("=")[1]
    MH = float(mass_over_z) * charge - (charge - 1) * 1.0078
    w_list = [title, charge_state, mass_over_z, str(MH)]
    frag_dic = {}
    for m in range(num1 + 5, num2):
        frag_mass = float(all[m].strip().split(" ")[0])
        frag_inten = all[m].strip().split(" ")[1]
        if frag_mass not in frag_dic:
            frag_dic[frag_mass] = float(frag_inten)
        else:
            print("waning: frag_dic wrong")
    frag_mass_list = sorted(list(frag_dic.keys()))
    frag_inten_list = list(frag_dic.values())

    inten_thd = max(frag_inten_list) * noise_cutoff
    scan_start = modi_mass_range[0] + feature_ion_mass_list[0]
    scan_end = modi_mass_range[1] + feature_ion_mass_list[0]
    for frag_mass in frag_mass_list:
        if frag_mass > min(frag_mass_list[-1], scan_end + 2) - 0.001:
            return False, ""
        else:
            if frag_mass <= scan_start - 2:
                continue
            elif frag_mass > scan_start - 2 and frag_mass < scan_end + 2:
                if frag_dic[frag_mass] > inten_thd:
                    feature_list = []
                    delta_mass = frag_mass - feature_ion_mass_list[0]
                    for word in feature_ion_mass_list:
                        feature_list.append(word + delta_mass)

                    judge_bool = judge_spectra(feature_list, frag_mass_list)
                    if judge_bool:
                        w_list.append(str(round(delta_mass, 2)))
                        return True, "\t".join(w_list)
                    else:
                        continue
                else:
                    continue


file_list = os.listdir(os.getcwd())
for name in file_list:
    if name[-3:] == "mgf":
        all = open(name).readlines()
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
            for i in range(len(spec_be_list)):
                judge_result, text_judege = get_spec_info(
                    spec_be_list[i], spec_end_list[i], all)

                if judge_result is False:
                    pass
                else:
                    b.write(text_judege)
                    b.write("\n")

        b.close()
    else:
        continue
