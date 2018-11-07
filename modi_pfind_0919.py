import os

# "KARGO_R_FREE[K]" "KARGO_R_H2O[K]" KARGO_K_FREE[R], KARGO_R_pH2O[R], KARGO_R_mH2O[R]
os.chdir(r"D:\workspace\pFindTask36_alpk1\result")

f = open("pFind.spectra", "r").readlines()
final_dic = {"ALPK1_NC": 0, "ALPK1_N": 0}

modifi = "KARGO_R_H2O[K]"
tgt_site = 383
n = 0
for line in f[1:]:
    if modifi in line:
        line_list = line.strip().split("\t")
        modi_list = line_list[10].split(";")[:-1]
        pep_start = line_list[13].split(",")[0]
        pep = line_list[5]
        for mo in modi_list:
            if modifi in mo:
                modi_on_pep = mo.split(",")[0]
                modi_on_pro = int(pep_start) + int(modi_on_pep)
                if modi_on_pro == tgt_site:                   
                    if pep[int(modi_on_pep)-1] == "K":
                        pass
                    else:
                        print(f.index(line))
                        print("wrong")
                    if "ALPK1_NC_" in line:
                        final_dic["ALPK1_NC"] += 1
                    elif "ALPK1_N_" in line:
                        final_dic["ALPK1_N"] += 1
                    else:
                        print(f.index(line)+"wrong")
                else:
                    continue
            else:
                continue

print(final_dic)
