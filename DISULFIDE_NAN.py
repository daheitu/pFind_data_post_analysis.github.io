import os

os.chdir(r"E:\CY\test")
f = open("N15_CROSS_LINK_SEPC.TXT", 'r').readlines()

pep_ratio_dic = {}
for line in f[1:]:
    line_list = line.strip().split("\t")
    pep = line_list[1]
    ratio = line_list[9]
    if pep not in pep_ratio_dic:
        pep_ratio_dic[pep] = [ratio]
    else:
        pep_ratio_dic[pep].append(ratio)

# print(pep_ratio_dic)
b = open("report.txt", 'w')

pep_rep_dic = {}
for pep in pep_ratio_dic:
    ratio_list = pep_ratio_dic[pep]
    NaN_ratio = ratio_list.count('9.77E-04')/len(ratio_list)
    pep_rep_dic[pep] = [pep, str(NaN_ratio), str(len(ratio_list))]
    b.write("\t".join(pep_rep_dic[pep]))
    b.write("\n")

b.close()