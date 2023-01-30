import os

f = open(r"E:\pFindWorkspace\TMT_test\result\pFind.protein").readlines()

pep_mod_dic = {}
# pep_mod_list = []
for line in f[2:]:
    linelist = line.split("\t") 
    if len(linelist) == 19:
        # print(linelist)
        pep = linelist[3]
        mods = linelist[8]
        spec_num = int(linelist[-1][:-1])
        # pep_mod_list.append((pep, mods, spec_num))
        pep_mod_dic[(pep, mods)] = spec_num


def get_TMT_ratio(pep_mod_dic):
    tmt_p = 0
    no_tmt_p = 0
    tmt_p_spec = 0
    no_tmt_p_spec = 0
    for pep, mod in pep_mod_dic:
        spec_num = pep_mod_dic[pep, mod]
        if "TMT6plex" not in mod:
            no_tmt_p += 1
            no_tmt_p_spec += spec_num
        else:
            tmt_p += 1
            tmt_p_spec += spec_num
    ratio_pep = tmt_p/(tmt_p + no_tmt_p)
    ratio_spec = tmt_p_spec/(tmt_p_spec + no_tmt_p_spec)
    return ratio_pep, ratio_spec


def get_total_TMTspec(tgt_pep, tmt_dic):
    tot_spec = 0
    for pep, mod in tmt_dic:
        if pep == tgt_pep:
            tot_spec += tmt_dic[pep, mod]
    
    return tot_spec


def get_tmt_sub_ratio(pep_mod_dic):
    b = open("no_tmt.csv", 'w')
    tmt_dic = {}
    no_tmt_dic = {}
    for pep, mod in pep_mod_dic:
        spec_num = pep_mod_dic[pep, mod]
        if "TMT6plex" not in mod:
            no_tmt_dic[pep, mod] = spec_num
        else:
            tmt_dic[pep, mod] = spec_num
    
    for pep, mod in no_tmt_dic:
        spec_no = no_tmt_dic[pep, mod]
        spec_tmt = get_total_TMTspec(pep, tmt_dic)
        print(pep, spec_no, spec_tmt)
        b.write("%s,%d,%d\n" % (pep, spec_no, spec_tmt))
    
    b.close()

print(get_TMT_ratio(pep_mod_dic))
print(get_tmt_sub_ratio(pep_mod_dic))
# for pep_info in pep_mod_list:
#     pep_mod = pep_info[:2]
#     spec = pep_info[-1]
#     if pep_mod not in pep_mod_dic:
#         pep_mod_dic[pep_mod] = [spec]
#     else:
#         pep_mod_dic[pep_mod].append(spec)

# for pep_mod in pep_mod_dic:
#     spec_list = pep_mod_dic[pep_mod]
#     if len(spec_list) > 1:
#         print(pep_mod, spec_list)

