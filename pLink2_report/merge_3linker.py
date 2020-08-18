# coding = utf-8

import os

wkdir = r"G:\msData\20200425\Lysozyme"


def get_link_sites(rpfl_path):
    sites_list = []
    f = open(rpfl_path).readlines()
    for line in f[1:]:
        lineList = line.split(",")
        if lineList[0].isdigit():
            break
        else:
            sites_list.append(lineList[0])
    return sites_list


rep_dic = {}
for root, dirs, fls in os.walk(wkdir):
    if root.endswith("reports"):
        linker, out_put = root.split("\\")[-3:-1]
        if linker == "DSSO":
            if "score" in root:
                for fl in fls:
                    if fl.endswith("v5.csv"):
                        rep_path = os.path.join(root, fl)
                        rep_dic[linker] =get_link_sites(rep_path)
        else:
            for fl in fls:
                if fl.endswith("v5.csv"):
                    rep_path = os.path.join(root, fl)
                    rep_dic[linker] =get_link_sites(rep_path)


b = open(wkdir.split("\\")[-1]+"_summary.csv", 'w')
print(rep_dic)
row = len(rep_dic["DSS"])
for i in range(row):
    # if i > len(rep_dic[])
    wlist = []
    for linker in ["DSS", "BSMEG", "DSSO"]:
        if i >= len(rep_dic[linker]):
            wlist.append("")
        else:
            wlist.append(rep_dic[linker][i])
    b.write(",".join(wlist)+"\n")
b.close()