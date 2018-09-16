import os
spec_cutoff = 2
sample_list = [
    "BSA_argo2_Tirs(BORATE)", "BSA_argo2_Tirs", "BSA_argo2_borae",
    "BSA_argo2_hepes(BORATE)", "BSA_argo2_hepes"
]

sample_list_GST = [
    "GST_50mM_BH", "GST_50mM_BTirs", "GST_50mM_Borate", "GST_50mM_HEPES",
    "GST_50mM_Tirs"
]
sample_list_aldo = [
    "Aldonase_argo2_borae", "Aldonasea_argo2_Tirs(BORATE)",
    "Aldonasea_argo2_Tirs", "Aldonasea_argo2_hepes(BORATE)",
    "Aldonasea_argo2_hepes"
]
os.chdir(
    r"C:\Users\Yong\Documents\pLink\pLink_task_2018.05.22.15.46.42_BSA_BUFFER_SCREEN\reports"
)
pep_pos_list = []
for fl in os.listdir(os.getcwd()):
    if "cross-linked_peptides" in fl:
        f = open(fl, 'r').readlines()
        for i in range(2, len(f)):
            line_list = f[i].strip().split(",")
            if len(line_list) == 6:
                pep_pos_list.append(i)
            else:
                continue
        pep_pos_list.append(len(f))
    else:
        continue

pep_dic = {}
for samp in sample_list:
    pep_dic[samp] = 0
    for i in range(len(pep_pos_list) - 1):
        low_pos = pep_pos_list[i]
        high_pos = pep_pos_list[i + 1]
        if high_pos - low_pos < spec_cutoff:
            continue
        else:
            for j in range(low_pos, high_pos):
                if samp in f[j]:
                    pep_dic[samp] += 1
                    break
                else:
                    continue
print(pep_dic)

b = open("peptide_rep_2.txt", 'w')
for word in pep_dic:
    w_list = [word, str(pep_dic[word])]
    b.write("\t".join(w_list))
    b.write("\n")
b.close()
