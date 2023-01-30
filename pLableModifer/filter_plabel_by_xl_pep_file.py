
import imp


import os

from scipy.fftpack import tilbert


fL_path = r"Z:\STY_PROJ\median_comp_sample\ribosome_Ecoli\pLink\Ribosome_for_XL_all_contained_Proteins_DSSO_GVL_output_plink_dsso\reports\filter_xlpep.csv"
plable_path = r"Z:\STY_PROJ\median_comp_sample\ribosome_Ecoli\pLink\Ribosome_for_XL_all_contained_Proteins_DSSO_GVL_output_plink_dsso"

raw_spec_dic = {}
f = open(fL_path).readlines()
for line in f[:]:
    linelist = line.split(",")
    if linelist[0] == "":
        title = linelist[2]
        raw_name = title.split(".")[0]
        # print(linelist[2].split(".")[0])
        if raw_name not in raw_spec_dic:
            raw_spec_dic[raw_name] = [title]
        else:
            raw_spec_dic[raw_name].append(title)

print(raw_spec_dic)


def in_spec_list(title, spec_list):
    for spec in spec_list: 
        if spec.upper() == title:
            return True
    return False

def get_palebel_path(plabel_root, rawname):
    for fl in os.listdir(plabel_root):
        if fl.startswith(rawname) and fl.endswith(".plabel") and "filter" not in fl and ".cross-linked" in fl:
            return os.path.join(plabel_root, fl)


for raw in raw_spec_dic:
    coor_plable_path = get_palebel_path(plable_path, raw)
    spec_list = raw_spec_dic[raw]
    print(spec_list)
    print(coor_plable_path)
    f = open(coor_plable_path).readlines()
    tgt_path = os.path.join(os.path.dirname(coor_plable_path), os.path.basename(coor_plable_path)[:-7]+"filter_spe.plabel")
    b = open(tgt_path, 'w')
    for i in range(8):
        b.write(f[i])
    wlist = []
    n = 0
    for i in range(9, len(f), 3):
        title = f[i+1].split("=")[1].strip()
        # print(title)
        in_s = in_spec_list(title, spec_list)
        if in_s:
            n += 1
            wlist.append("[Spectrum%d]\n" % n)
            wlist.append(f[i+1])
            wlist.append(f[i+2])
    wlist.insert(0, "total=%d\n" % n)
    b.writelines(wlist)
    b.close()



