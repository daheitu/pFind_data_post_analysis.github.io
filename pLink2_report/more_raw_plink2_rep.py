# coding = utf-8

from plink2_report_v5 import main_flow
import os
import re


wd_dir = r"G:\msData\20200419\BSA"

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
                prot, linker, out = dir_list[-4:-1]
                if fl.endswith("v5.csv") and "pep" not in fl:
                    sample_name = "_".join([prot, linker, out])
                    f = open(os.path.join(root, fl)).readlines()
                    for line in f:
                        linelist = line.split(",")
                        if linelist[0].isdigit():
                            site_num, spec_num = linelist[:2]
                    data_dic[sample_name] = [sample_name, site_num, spec_num]

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
        prot, linker = name_process(key)
        if "mass_filter" not in linker:
            wlist = [prot, linker]
            wlist.extend(dic[key][-2:])
            b.write(",".join([str(x) for x in wlist])+"\n")
    b.close()

if __name__ == "__main__":
    summary()
    rep_dic = merge_data()
    writedic2file(rep_dic, "rep.csv")
    print(os.getcwd())