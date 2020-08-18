path1 = r"C:\Users\Yong Cao\Documents\pLink\LuShan_NP\reports\LuS_NP_2020.08.02_BS3_Trypsin_v5.csv"
path2 = r"C:\Users\Yong Cao\Documents\pLink\Lushan_np_bs3_d4\reports\LuS_NP_2020.08.02_BS3_heavy_Trypsin_v5.csv"

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
            info_dic[site] = [spec, evlue, svm]
    return info_dic


info_d0_dic = get_info(path1)
info_d4_dic = get_info(path2)
site_d0_set = set(list(info_d0_dic.keys()))
site_d4_set = set(list(info_d4_dic.keys()))

b = open("report.csv", 'w')
for site in list(site_d0_set| site_d4_set):
    if site in info_d0_dic and site in info_d4_dic:
        wlist = [site]
        wlist.extend(info_d0_dic[site])
        wlist.extend(info_d4_dic[site])
    elif site in info_d0_dic:
        wlist = [site]
        wlist.extend(info_d0_dic[site])
    else:
        wlist = [site]
        wlist.extend(["", "", ""])
        wlist.extend(info_d4_dic[site])
    b.write(",".join(wlist)+"\n")
b.close()

