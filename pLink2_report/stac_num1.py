# coding = utf-8
import os
wd_dir = r"K:\LiuFan_FAIMS_AC\DSS\FAIMS_DSS_LIU\reports"
def stac_num1(fl_path):
    site1 = 0
    f = open(fl_path).readlines()
    for line in f[1:]:
        linelist = line.split(",")
        if linelist[0].isdigit():
            tot_sites = int(linelist[0])
            break
        else:
            spec = linelist[1]
            if spec == "1":
                site1 += 1

    return site1/tot_sites


for root, dirs, fls in os.walk(wd_dir):
    dir_list = root.split("\\")
    if dir_list[-1] == "reports":
        print(root)
        for fl in fls:
            prot, linker, out_type = dir_list[-4:-1]
            if fl.endswith("v5.csv") and "pep" not in fl:
                num1_ratio = stac_num1(os.path.join(root, fl))
                print(num1_ratio)