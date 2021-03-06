import os
os.chdir(
    r"F:\ALL ArGO data\ArGO\optimize condition\Buffer screen\BSA\blank_ArGO")
spectra_cutoff = 4


def get_info(opened_file):
    modi_dic = {}
    for line in opened_file[3:]:
        if line[:3] == "---":
            break
        else:
            line_list = line.strip().split("\t")
            if len(line_list) == 25:
                peptide = line_list[1]
                modificaiton = line_list[6]
                spec_num = int(line_list[-1])
                if "PFIND_DELTA" in modificaiton:
                    modi_list = modificaiton.split(";")[:-1]
                    for modi in modi_list:
                        if "PFIND_DELTA" in modi:
                            modi_pos = modi.split(",")[0]
                            delta_mass = modi.split(",")[1].split("_")[-1]
                            modi_site = peptide[int(modi_pos) - 1]
                            if modi_site not in modi_dic:
                                modi_dic[modi_site] = {}
                                if delta_mass not in modi_dic[modi_site]:
                                    modi_dic[modi_site][delta_mass] = spec_num
                                else:
                                    pass
                            else:
                                if delta_mass not in modi_dic[modi_site]:
                                    modi_dic[modi_site][delta_mass] = spec_num
                                else:
                                    modi_dic[modi_site][delta_mass] += spec_num
                        else:
                            continue
                else:
                    continue
            else:
                continue

    print(modi_dic)
    return modi_dic


def sort_str_digtal_list(list):
    final_list = []
    new_list = []
    for char in list:
        if char.isdigit():
            new_list.append(int(char))
        else:
            new_list.append(float(char))

    new_list.sort()
    for num in new_list:
        if len(str(num).split(".")[1]) == 1:
            final_list.append(str(num) + "0")
        else:
            final_list.append(str(num))

    return final_list


def main():
    f = open("pFind.protein", 'r').readlines()
    report_dic = get_info(f)
    b = open("report.txt", 'w')
    AAs = list(report_dic.keys())
    AAs.sort()
    for AA in AAs:
        b.write(AA)
        b.write("\n")
        mass_1 = list(report_dic[AA].keys())
        masses = sort_str_digtal_list(mass_1)
        for mass in masses:
            if report_dic[AA][mass] > spectra_cutoff:
                b.write("\t".join([mass, str(report_dic[AA][mass])]))
                b.write("\n")
            else:
                continue
    b.close()


if __name__ == "__main__":
    main()
