import re, os
wd_dir = r"G:\msData\SLu\G1"
# f = open(r"G:\msData\SLu\G1\LuS_NP_output_BS3\reports\pQuant_xl\pQuant.proteins.list").readlines()

titleLine = "Site,d0_spec,d0_best_evalue, d0_svm,d0_pep,d4_spec,d4_best_evalue, d4_svm,d4_pep, d4/d0 ratio, d4/d0 score\n"

def judgeHomoHetro(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]
    if int(position1) > int(position2):
        return protein2 + "(" + position2 + ")-" +  protein1 + "(" + position1 + ")"
    else:
        return linked_site



def get_quant_dic(wd_dir):
    site_dic = {}
    qt_path = wd_dir + "\\LuS_NP_output_BS3\\reports\\pQuant_xl\\pQuant.proteins.list"
    f = open(qt_path).readlines()
    for line in f[1:]:
        linelist = line.split("\t")
        site = linelist[0]
        d40_ratio = linelist[3]
        interf_score = linelist[4]
        site = judgeHomoHetro(site)
        if site not in site_dic:
            site_dic[site] = [[d40_ratio], [interf_score]]
        else:
            site_dic[site][0].append(d40_ratio)
            site_dic[site][1].append(interf_score)

    for site in site_dic:
        ratio = ";".join(site_dic[site][0])
        score = ";".join(site_dic[site][1])
        site_dic[site]= [ratio, score]


    return site_dic


def get_rp_fl_path(rep_dir):
    for fl in os.listdir(rep_dir):
        if fl.endswith("v5.csv"):
            return os.path.join(rep_dir, fl)


def get_rep_info_dic(fl_path):
    rep_dic = {}
    f = open(fl_path).readlines()
    for line in f[1:]:
        linelist = line.split(',')
        site = linelist[0]
        if site.isdigit():
            break
        else:
            rep_dic[site] = linelist[1:5]
    return rep_dic


def merge_d0d4(root_path):
    qt_dic = get_quant_dic(root_path)
    repd0_dir = root_path + "\\LuS_NP_output_BS3\\reports"
    repd4_dir = root_path + "\\LuS_NP_output_BS3_heavy\\reports"
    d0_path = get_rp_fl_path(repd0_dir)
    d4_path = get_rp_fl_path(repd4_dir)
    d0_info_dic = get_rep_info_dic(d0_path)
    d4_info_dic = get_rep_info_dic(d4_path)
    site_set = set(list(d0_info_dic.keys())) | set(list(d4_info_dic.keys()))
    b = open(os.path.join(root_path, "group%s_xlink.csv" % wd_dir[-1]), 'w')
    b.write(titleLine)
    for site in list(site_set):
        wlist = [site]
        if site in d0_info_dic and site in d4_info_dic:
            wlist.extend(d0_info_dic[site])
            wlist.extend(d4_info_dic[site])
            wlist.extend(qt_dic[site])
        elif site in d0_info_dic:
            wlist.extend(d0_info_dic[site])
            wlist.extend([""]*4)
            wlist.extend(qt_dic[site])
        else:
            wlist.extend([""]*4)
            wlist.extend(d4_info_dic[site])
        b.write(",".join(wlist)+"\n")
    b.close()

        
merge_d0d4(wd_dir)