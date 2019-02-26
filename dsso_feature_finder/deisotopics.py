import os

isotopic_tol = 20
os.chdir(r"/Users/yong/AnacondaProjects")
f = open("./deiso_spec.mgf", 'r').readlines()

mass_list = []
mass_tuple = tuple(mass_list)
ins_list = []
ms_ins_dic = {}
for line in f:
    if line[:5] == "TITLE":
        title = line[6:-1]
    elif line[0].isdigt:
        ms = line.strip().split(" ")[0]
        ins = line.strip().split(" ")[1]
        ms_ins_dic[ms] = ins
        mass_list.append(ms)
        ins_list.append(float(ins))
    elif line.strip() == "END IONS":
        break
    else:
        continue


sorted_int_list = sorted(ins_list, reverse=True)

mz_chrg = {}
for i range(len(sorted_int_list)):
    pos = ins_list.index(sorted_int_list[i])
    mz = mass_tuple[pos]
    if mz in mass_list:
        mz_next = mass_tuple[pos+1]
        if mz_next - mz > 1:
            mz_chrg[mz] = 0
        else:






    