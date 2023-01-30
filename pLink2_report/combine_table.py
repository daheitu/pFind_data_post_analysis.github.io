import os
path1 = r"M:\collibrations\YangSS_rpa_20211005\RPA_105_100NT\RTT105_RPA_con_EDC-DE_output\reports\RTT105_RPA_con_2021-10-05_EDC-DE_Trypsin_v5.csv" #结果文件1的路径
path2 = r"M:\collibrations\YangSS_rpa_20211005\RPA_100NT\RTT105_RPA_con_EDC-DE_output\reports\RTT105_RPA_con_2021-10-05_EDC-DE_Trypsin_v5.csv" # 结果文件2的路径

des_dir = r"D:\collaboration\yangSS_LIQlab\20210821"
file_name = "YangSS_RPAvsRPA_105_100NT_0426_cutoff.csv"
need_nromal= True

delta_spec_cutoff = 6


####################Don't change the following lines#########################
def get_info(path1):
    info_dic = {}
    f = open(path1).readlines()
    for line in f[1:]:
        linelist = line.split(",")
        if linelist[0].isdigit():
            break
        else:
            site = linelist[0]
            spec = linelist[1]
            evlue = linelist[2]
            svm = linelist[3]
            link_type = linelist[5]
            info_dic[site] = [spec, evlue, svm, link_type]
    return info_dic


def get_normal_ratio():
    if not need_nromal:
        return 1.00
    else:
        info_d0_dic = get_info(path1)
        info_d4_dic = get_info(path2)
        total_spec_f1 = sum([float(x[0]) for x in info_d0_dic.values()])
        total_spec_f2 = sum([float(x[0]) for x in info_d4_dic.values()])
        normal_ratio = round(total_spec_f1/total_spec_f2, 5)
        return normal_ratio 





def is_significant(spec1 ,spec2, normal_ratio, cutoff= 5):
    if "" in [spec1, spec2]:
        if spec1 == "":
            spec2 = float(spec2) * normal_ratio
        else:
            spec1 = float(spec1) / normal_ratio

        other = float([x for x in [spec1, spec2] if x != ""][0])
        if other < cutoff:
            return "no"
        else:
            if spec1 == "":
                return "up"
            else:
                return "down"
    else:
        spec2 = float(spec2) * normal_ratio
        min_one = min([float(x) for x in [spec1, spec2]])
        max_one = max([float(x) for x in [spec1, spec2]])
        if max_one <= float(cutoff) * 2:
            if max_one - min_one < cutoff:
                return "no"
            else:
                if float(spec1) > float(spec2):
                    return "down"
                else:
                    return "up"
        else:
            if max_one / min_one < 2:
                return "no"
            else:
                if float(spec1) > float(spec2):
                    return "down"
                else:
                    return "up"


def main():
    name1 = os.path.basename(path1)
    name2 = os.path.basename(path2)
    b = open(os.path.join(des_dir, file_name), 'w')
    b.write("site,%s,,,%s,,,inter or intra,is_significant\n"%(name1, name2))
    info_d0_dic = get_info(path1)
    info_d4_dic = get_info(path2)
    site_d0_set = set(list(info_d0_dic.keys()))
    site_d4_set = set(list(info_d4_dic.keys()))

    normal_ratio = get_normal_ratio()

    for site in list(site_d0_set| site_d4_set):
        wlist = [site]

        if site in info_d0_dic and site in info_d4_dic:
            wlist.extend(info_d0_dic[site][:-1])
            wlist.extend(info_d4_dic[site][:-1])
            wlist.append(info_d0_dic[site][-1])
            spec1 = info_d0_dic[site][0]
            spec2 = info_d4_dic[site][0]
        elif site in info_d0_dic:
            # wlist = [site]
            wlist.extend(info_d0_dic[site][:-1])
            wlist.extend(["", "", ""])
            wlist.append(info_d0_dic[site][-1])
            spec1 = info_d0_dic[site][0]
            spec2 = ""

        else:
            # wlist = [site]
            wlist.extend(["", "", ""])
            wlist.extend(info_d4_dic[site])
            spec1 = ""
            spec2 = info_d4_dic[site][0]
        
        wlist.append(is_significant(spec1, spec2, normal_ratio, delta_spec_cutoff))
        b.write(",".join(wlist)+"\n")
    b.write("The normalize ratio is %f\n" % normal_ratio)
    b.close()


if __name__ == "__main__":
    main()
    print("Well Done")