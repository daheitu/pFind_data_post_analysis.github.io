import os
import re
import statistics

MB_dif_cuoff = 0.15


def format_ratio(numb):
    if numb < 0.001:
        return str(format(numb, '.2e'))
    else:
        return str(round(numb, 3))


def get_char_pos_list(str2, str1):
    pos_list = []
    if str1.count(str2) == 0:
        return pos_list
    else:
        for i in range(len(str1)):
            if str1[i] == str2:
                pos_list.append(i)
            else:
                continue
        return pos_list


def get_linked_site_inform(linked_site):
    print(linked_site)
    pos_list = re.findall("\((\d*)\)", linked_site)
    return pos_list


def judgeAndCal(link_site, ratio_list_mono, ratio_list_best):
    pos1 = get_linked_site_inform(link_site)[0]
    pos2 = get_linked_site_inform(link_site)[1]
    mono_list = ratio_list_mono.split(";")
    best_list = ratio_list_best.split(";")
    noma_list = []
    for i in range(len(mono_list)):
        noma_list.append((float(mono_list[i]) + float(best_list[i])) / 2)
    Inter_pect = (noma_list[1] + noma_list[0]) / (1 + noma_list[2])
    if Inter_pect == 1:
        return "Inter", 1000.0
    else:
        Ratio_Inter_Intra = Inter_pect / (1 - Inter_pect)
        if abs(int(pos1) - int(pos2)) < 5:
            return "Inter", Ratio_Inter_Intra
        else:
            if Ratio_Inter_Intra < 0:
                return "Inter", Ratio_Inter_Intra
            elif Ratio_Inter_Intra < 0.333:
                return "Intra", Ratio_Inter_Intra
            elif Ratio_Inter_Intra > 3:
                return "Inter", Ratio_Inter_Intra
            else:
                return "Mixture", Ratio_Inter_Intra


def simplify_protein_list(f):
    site_dic = {}
    n15_ratio_list = []
    for line in f:
        if line.rstrip("\n") != "":
            line_list = line.strip().split("\t")
            if ";" in line_list[-2]:
                n15 = line_list[-2][:-1].split(';')
                for ratio_sigma in n15:
                    n15_ratio_list.append(float(ratio_sigma.split(',')[0]))
        else:
            print(line)
    if n15_ratio_list == []:
        return {}
    else:
        mdi = statistics.median(n15_ratio_list)

        for line in f:
            line_list = line.strip().split("\t")
            if line_list[3][0].isdigit():
                ratio_ori = ";".join(
                    [line_list[3], line_list[10], line_list[17]])
                ratio_cor = ";".join([
                    format_ratio(float(line_list[3]) / mdi),
                    format_ratio(float(line_list[10]) / mdi),
                    format_ratio(float(line_list[17]) / mdi)
                ])
                sigma = ";".join([line_list[4], line_list[11], line_list[18]])
                site_dic[line_list[0]] = [
                    line_list[0], line_list[2], line_list[8], ratio_ori,
                    ratio_cor, sigma
                ]
        return site_dic


def cydiv(str1, str2):
    num1 = float(str1)
    num2 = float(str2)
    if num1 != 0.0 and num2 != 0.0:
        if num1 > 0.2 and num2 > 0.2:
            if abs(num1 - num2) / max(num1, num2) < MB_dif_cuoff:
                return True
            else:
                return False
        else:
            if abs(num1 - num2) < MB_dif_cuoff - 0.05:
                return True
            else:
                return False
    elif num1 == 0.0 and num2 == 0.0:
        return True
    else:
        if abs(num1 - num2) < MB_dif_cuoff - 0.05:
            return True
        else:
            return False


"D:\E\Collabaration\TC\DIUB_5\Quant\B\EGS\quant\Mono"
"pQuant.proteins_EGS.list"

# for char in ["A"]:  #"B" , "C", "D",
for linker in ["BS3", "BS2G", "DSS"]: # , "EGS"]
    path = "D:\\E\\Collabaration\\TC\DiUB_RNP\\" + "\\" + linker + "\\" + "quant"
    print(path)
    os.chdir(path)
    report_file = "DiUb"  + "_" + linker + "_" + "QUANT.txt"
    b = open(report_file, 'w')
    b.write("\t".join([
        "Link_site", "Leading group", "PSM", "Mono ratio origin",
        "Mono ratio Correction", "Mono sigma", "Best ratio origin",
        "Best ratio Correction", "Best sigam"
    ]))
    b.write("\n")

    file1_name = path + "\\Mono\\pQuant.proteins_" + linker + ".list"
    file2_name = path + "\\Best\\pQuant.proteins_" + linker + ".list"
    f1 = open(file1_name, 'r').readlines()
    f2 = open(file2_name, 'r').readlines()
    print(len(f1), len(f2))
    f1_dic = simplify_protein_list(f1)
    f2_dic = simplify_protein_list(f2)
    if f1_dic == {} and f2_dic == {}:
        continue
    else:
        for site in f1_dic:
            if site in f2_dic.keys():
                if f1_dic[site][2] == f2_dic[site][2] and f1_dic[site][1] == f2_dic[site][1]:
                    f1_dic[site].extend(f2_dic[site][3:6])

                else:
                    print("spectra numer wrong" + site)
            else:
                print("file wrong")
            
            mono_ratio_list = f1_dic[site][4].split(";")
            Best_ratio_list = f1_dic[site][7].split(";")

            AlBh = cydiv(mono_ratio_list[0], Best_ratio_list[0])
            AhBl = cydiv(mono_ratio_list[1], Best_ratio_list[1])
            AhBh = cydiv(mono_ratio_list[2], Best_ratio_list[2])
            
            if [AlBh, AhBl, AhBh].count(True) == 3:
                f1_dic[site].append("right")
                Link_type, Inter_Intra_ratio = judgeAndCal(
                    f1_dic[site][0], f1_dic[site][4], f1_dic[site][7])
                f1_dic[site].append(",".join([
                    Link_type,
                    "(" + str(format_ratio(Inter_Intra_ratio)) + ")"
                ]))
            else:
                f1_dic[site].append("wrong")
                f1_dic[site].append("")

            f1_dic[site]
            b.write("\t".join(f1_dic[site]))
            b.write("\n")
    print("done")
    b.close()

print("Done")
