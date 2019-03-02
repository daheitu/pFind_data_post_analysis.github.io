pep = 'D(V/L/A)KPSTEHI(H/W)LK'

import re
import copy
pep = 'D(V/L/A)KPSTEHI(H/W)LK'
f = open('./peptides_sythensis_sgc.csv', 'r').readlines()
resi_num_list = [38]




AA_resi_dic = {'A':71.037114,'R':156.101111,'N':114.042927,'D':115.026943,'C':103.009185, \
'E':129.042593,'Q':128.058578,'G':57.021464,'H':137.058912,'I':113.084064, \
'L':113.084064,'K':128.094963,'M':131.040485,'F':147.068414,'P':97.052764, \
'S':87.032028,'T':101.047679,'U':150.95363,'W':186.079313,'Y':163.06332,'V':99.068414, \
'H2O':18.01056,'Proton':1.0072766}


def cal_pep_mass(seq):
    mass = AA_resi_dic["H2O"]
    for resi in list(seq):
        mass += AA_resi_dic[resi]
    return round(mass, 3)

def deal_vari_pep(pep):
    a = re.findall("\((.*?)\)", pep)
    b = re.split("\((.*?)\)", pep)
    # print(b)
    rep_dic = {}
    a_len = len(a)
    if a_len == 0:
        print("wrong")
    elif a_len == 1:
        vair_1_list = a[0].split('/')
        vair_idx = b.index(a[0])
        new_list = b.copy()
        for char in vair_1_list:
            new_list[vair_idx] = char
            new_seq = "".join(new_list)
            new_seq_mass = cal_pep_mass(new_seq)
            rep_dic[new_seq] = new_seq_mass

    elif a_len == 2:
        vair_1_list = a[0].split('/')
        vair_2_list = a[1].split('/')
        vair1_idx = b.index(a[0])
        vair2_idx = b.index(a[1])       
        for char in vair_1_list:
            for char2 in vair_2_list:
                new_list = b.copy()
                new_list[vair1_idx] = char
                new_list[vair2_idx] = char2
                new_seq = "".join(new_list)
                new_seq_mass = cal_pep_mass(new_seq)
                rep_dic[new_seq] = new_seq_mass
    else:
        print("please contact Yong Cao")
    return rep_dic

# print(rep_dic)
aq_pep_list = []
resi_pep_list = []
for line in f[1:]:
    line_list = line.strip().split(',')
    pep = line_list[4]
    name = line_list[1]
    name_num = int(name[2:])
    if name_num not in resi_num_list:
        aq_pep_list.append([name, pep])
    else:
        pass
        # resi_pep_list.append([name, pep])

print(aq_pep_list)

resi_pep_list = aq_pep_list.copy()
print(resi_pep_list)
final_dic = {}
for aq_pep in aq_pep_list:
    name_aq = aq_pep[0]
    pep_aq = aq_pep[1]
    pep_aq_dic = deal_vari_pep(pep_aq)
    for resi_pep in resi_pep_list:
        name_re = resi_pep[0]
        pep_re = resi_pep[1]
        pep_re_dic = deal_vari_pep(pep_re)
        combin_dic = {}
        for sPep_aq in pep_aq_dic:
            for sPep_re in pep_re_dic:
                combin_pep = (sPep_aq, sPep_re)
                combin_mass = pep_aq_dic[sPep_aq] + pep_re_dic[sPep_re]
                combin_dic[combin_pep] = combin_mass
        mass_list = list(combin_dic.values())
        mass_list.sort()
        delta_list = []
        for i in range(1, len(mass_list)):
            delta_list.append(mass_list[i] - mass_list[i -1])
        dtm_m4 = 0
        for dtm in delta_list:
            if dtm < 4:
                dtm_m4 += 1
        ave_delta_mass = sum(delta_list)/len(delta_list)
        final_dic[(name_aq, name_re)] = [combin_dic, ave_delta_mass, \
            min(delta_list), max(delta_list), dtm_m4, len(delta_list)+1]

print(final_dic)
b = open('report_test_all2all.csv', 'w')
for pair in final_dic:
    b.write(",".join([str(pair), str(final_dic[pair][1]), \
        str(final_dic[pair][2]), str(final_dic[pair][3]), \
        str(final_dic[pair][4]), str(final_dic[pair][5])]) + "\n")
    for pep_pair in final_dic[pair][0]:
        b.write(",".join(["", str(pep_pair), str(final_dic[pair][0][pep_pair])] )+ "\n")
b.close()
