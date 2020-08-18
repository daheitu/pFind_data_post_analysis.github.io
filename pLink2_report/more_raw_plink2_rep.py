# coding = utf-8

from plink2_report_v5 import main_flow
import os
import re


wd_dir = r"K:\20200801"

spec_cutoff = 3  # spectra number cut-off
Best_evalue_cutoff = 2 # 交联位点对层次最好的e-value cutoff
E_value_cutoff_SpecLvl = 2 # 谱图层次的e-value cutoff


def summary():
    for root, dirs, fls in os.walk(wd_dir):
        # print(root)
        if root.split("\\")[-1] == "reports":
            print(root)
            main_flow(root, spec_cutoff, Best_evalue_cutoff)


def merge_data():
    data_dic = {}
    for root, dirs, fls in os.walk(wd_dir):
        dir_list = root.split("\\")
        if dir_list[-1] == "reports":
            print(root)
            for fl in fls:
                prot, linker, out_type = dir_list[-4:-1]
                if fl.endswith("v5.csv") and "pep" not in fl:
                    print(prot, linker)
                    sample_name =  "_".join([prot, linker, out_type])
                    f = open(os.path.join(root, fl)).readlines()
                    n_intra = 0
                    n_inter = 0
                    n_other = 0
                    for line in f[1:]:
                        linelist = line.split(",")
                        if linelist[0].isdigit():
                            site_num, spec_num = linelist[:2]
                            break
                        else:
                            pro_type = linelist[5]
                            if pro_type == "Intra":
                                n_intra += 1
                            elif pro_type == "Inter":
                                n_inter += 1
                            else:
                                n_other += 1
                        
                    data_dic[sample_name] = [prot, linker, out_type, site_num, spec_num, n_intra, n_inter, n_other]

    return data_dic


def name_process(key_name):
    nlist = key_name.split('_')
    prot, linker = nlist[:2]
    if linker == "DSSO":
        stype = "_".join(nlist[-2:])
        linker = "_".join([linker, stype])
    return prot, linker


def writedic2file(dic, file_name):
    b = open(file_name, 'w')
    for key in dic:
        # linker, type_score = key.split('_output')
        # prot, linker = name_process(key)
        # if "mass_filter" not in key:
        wlist = dic[key]
        # wlist.extend(dic[key][-2:])
        b.write(",".join([str(x) for x in wlist])+"\n")
    b.close()


if __name__ == "__main__":
    summary()
    rep_dic = merge_data()
    print(rep_dic)
    writedic2file(rep_dic, "yeast_ribosome_sc_cf%d.csv" % spec_cutoff)
    print(os.getcwd())