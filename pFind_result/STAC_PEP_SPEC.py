import os
wd = r"M:\LiuXT\DSS_PFIND\pFindTask8\result"

modi = "Xlink_DSS[156][K]"

def get_protein_fl(wd):
    for fl in os.listdir(wd):
        if fl.endswith(".protein"):
            return os.path.join(wd, fl)
    return None

def is_decoy_pep(pro):
    pro_list = pro[:-1].split("/")
    for pro in pro_list:
        if "REV_" not in pro:
            return False
    return True


print(get_protein_fl(wd))
pfind_path = get_protein_fl(wd)
total_spec = 0
total_pep = 0

f = open(pfind_path).readlines()
for line in f[2:]:
    linelist = line.split("\t")
    if len(linelist) == 19 and modi in line:
        spec = int(linelist[-1][:-1])
        pros = linelist[10]
        if not is_decoy_pep(pros):
            total_spec += spec
            total_pep += 1

print("%d\t%d" % (total_pep, total_spec))
