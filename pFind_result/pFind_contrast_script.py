# coding = utf-8
import copy, os


wd_path = r"E:\workspace\pFindTask97_E.coli\result" # important the pFind.protein path
output_name = "Modification_filter.txt" # important! the output file name
modify = "Carbamidomethyl[C]" # target modification name


def get_modi_info(linelist, mod_tgt):
    pros = linelist[10]
    mods = linelist[8]
    pep = linelist[3]
    pros_list = pros[:-1].split('/')
    pep_pos_list = linelist[11][:-1].split('/')
    mod_list = mods[:-1].split(";")
    spec_num = linelist[-1]
    score = linelist[7]
    # spec_title = linelist[-3]
    info_list = []
    for mod_info in mod_list:
        md_pos, md_name = mod_info.split(',')
        if md_name  == mod_tgt:
            pro_pos_list = []
            for i in range(len(pros_list)):
                pro = pros_list[i]
                pep_pos = int(pep_pos_list[i].split(',')[0])
                mod_pos = pep_pos + int(md_pos)
                pro_pos = pro + "[" + str(mod_pos) + "]"
                pro_pos_list.append(pro_pos)
            pro_poss = "/".join(pro_pos_list) + "/"
            info_list.append([pros, pro_poss, pep, mods, score, spec_num])
    
    return info_list


#求互补蛋白序列                
def get_complement_proteins(protein_list, pro, group_pro):
    if len(protein_list) == 1 and protein_list[0] == pro:
        return ""
    else:
        group_set = set(group_pro) # group proteins
        pep_related_pro = set(protein_list)
        target_pro_set = set(pro)
        complement_pros = pep_related_pro - target_pro_set & group_set
        return "/".join(list(complement_pros))

    
def get_info_dic(fl_path):
    rep_dic = {}

    f = open(fl_path).readlines()
    i = 0
    while i < len(f):
        linelist = f[i].rstrip("\n").split("\t")
        if not linelist[0].isdigit():
            i += 1
        else:
            group_pro = [linelist[1]]
            # score_pro = linelist[]
            i += 1
            while i < len(f):
                sub_list = f[i].rstrip("\n").split("\t")
                if sub_list[2].isdigit():
                    break
                else:
                    group_pro.append(sub_list[2])
                    i += 1
            while i < len(f):
                pep_info_list = f[i].rstrip("\n").split("\t")
                if pep_info_list[0].isdigit():
                    break
                else:
                    pep =pep_info_list[3]
                    mod = pep_info_list[8]
                    pros = pep_info_list[10]
                    spec = pep_info_list[-1]
                    pro_list = pros[:-1].split("/")
                    best_score = pep_info_list[7]
                    is_unique = "not unique"
                    if len(pro_list) == 1:
                        is_unique = "is unique"
                    for pro in pro_list:
                        if pro in group_pro:
                            comple_pros = get_complement_proteins(pro_list, pro, group_pro)
                            if pro not in rep_dic:
                                rep_dic[pro] = {}
                                rep_dic[pro][pep, mod] = [spec, best_score, is_unique, comple_pros]
                            else:
                                rep_dic[pro][pep, mod] = [spec, best_score, is_unique, comple_pros]

                    i += 1
    return rep_dic


def reformate_dic(uni_dic):
    re_uni_dic = {}
    for pro in uni_dic:
        pro_p_num = 0
        pro_pep_unique = 0
        pro_spec_num = 0
        pro_spec_unique = 0
        poss_pep_dic = uni_dic[pro]
        have_uni_pep_pro = "no unique peptide"
        for pep_info in poss_pep_dic:
            spec, best_score, is_unique, comple_pros = poss_pep_dic[pep_info]
            pro_p_num += 1
            pro_spec_num += spec
            if is_unique == "is unique":
                pro_pep_unique += 1
                pro_spec_unique += spec
                have_uni_pep_pro= "have unique peptide"
        re_uni_dic[pro] = [["%d(%d)" % (pro_p_num, pro_pep_unique), "%d(%d)" % (pro_spec_num, pro_spec_unique),have_uni_pep_pro], poss_pep_dic]
        
        
    return re_uni_dic


def get_all_pros_sites_peps(total_dic):
    struced_dic = {}
    for sp in total_dic:
        re_uni_dic = total_dic[sp]
        for pro in re_uni_dic:
            if pro not in struced_dic:
                struced_dic[pro] = [[], {}]
            pro_sites_dic = re_uni_dic[pro][1]
            
            for site in pro_sites_dic:
                if site not in struced_dic[pro][1]:
                    struced_dic[pro][1][site] = [[], {}]
        
                site_peps_dic = pro_sites_dic[site][1]
                
                for pep_m in site_peps_dic:
                    struced_dic[pro][1][site][1][pep_m] = []
    return struced_dic


def compare_results(total_dic):
    if len(total_dic) == 1:
        return list(total_dic.values())[0]
    else:
        structed_dic = get_all_pros_sites_peps(total_dic)
        sp_list = list(total_dic.keys())
        for pro in structed_dic:
            for sp in sp_list:
                uni_dic = total_dic[sp]
                if pro not in uni_dic:
                    structed_dic[pro][0].extend(["", "", "", ""])
                    for site in structed_dic[pro][1]:
                        structed_dic[pro][1][site][0].extend(["", "", "", ""])

                        for pep_m in structed_dic[pro][1][site][1]:
                            structed_dic[pro][1][site][1][pep_m].extend(["", "", "", ""])
                else:
                    structed_dic[pro][0].extend(uni_dic[pro][0])
                    for site in structed_dic[pro][1]:
                        if site not in uni_dic[pro][1]:
                            structed_dic[pro][1][site][0].extend(["", "", "", ""])
                            
                            for pep_m in structed_dic[pro][1][site][1]:
                                structed_dic[pro][1][site][1][pep_m].extend(["", "", "", ""])

                        else:
                            structed_dic[pro][1][site][0].extend(uni_dic[pro][1][site][0])
                            for pep_m in structed_dic[pro][1][site][1]:
                                if pep_m not in uni_dic[pro][1][site][1]:
                                    structed_dic[pro][1][site][1][pep_m].extend(["", "", "", ""])
                                else:
                                    structed_dic[pro][1][site][1][pep_m].extend(uni_dic[pro][1][site][1][pep_m])
        return structed_dic


def write_dic2_file(final_dic, sp_list, w_name):
    b = open(os.path.join(wd_path, w_name), 'w')
    line1_list = [""] * 4
    for sp in sp_list:
        line1_list.extend([sp, "", "", ""])
    b.write("\t".join(line1_list)+"\n")
    line2_list = ["Protein", "", "", ""]
    line2_list.extend(["Total_site_num@pro","Total_pep_num@pro", "Total_spec_num@pro", "have unique peptide?"]* len(sp_list))
    b.write("\t".join(line2_list)+"\n")
    line3_list = ["", "site", "", ""]
    line3_list.extend(["Total_pep_num@site(Total_unique_pep_num@site)", "Total_spec_num@site(Total_unique_spec_num@site)", "best-score@site", "have unique peptide?"]* len(sp_list))
    b.write("\t".join(line3_list)+"\n")
    line4_list = ["", "", "Peptide", "Modification"]
    line4_list.extend(["Total_spec_num@pep", "shared proteins", "best-score@pep", "is unique?"]* len(sp_list))
    b.write("\t".join(line4_list)+"\n")

    for pro in final_dic:
        pro_info, sites_dic = final_dic[pro]
        pro_wlist = [pro, "", "", ""]
        pro_wlist.extend(pro_info)
        b.write("\t".join([str(x) for x in pro_wlist]) + "\n")
        for site in sites_dic:
            site_info, peps_dic =  sites_dic[site]
            site_wlist = ["", site, "", ""]
            site_wlist.extend(site_info)
            b.write("\t".join([str(x) for x in site_wlist]) + "\n")
            for pep, mod in peps_dic:
                pep_wlist = ["", "", pep, mod]
                pep_wlist.extend(peps_dic[pep,mod])
                b.write("\t".join([str(x) for x in pep_wlist]) + "\n")
    b.close()



def main():
    total_dic = {}
    for fl in os.listdir(wd_path):
        if fl.endswith(".protein"):
            fl_path = os.path.join(wd_path, fl)
            fl_name = fl[:-8]
            total_dic[fl_name] =  reformate_dic(get_info_dic(fl_path))
    final_dic = compare_results(total_dic)
    sp_list = list(total_dic.keys())
    print(final_dic)
    write_dic2_file(final_dic, sp_list, output_name)
    


if __name__ == "__main__":
    main()
    print("Well Done!")