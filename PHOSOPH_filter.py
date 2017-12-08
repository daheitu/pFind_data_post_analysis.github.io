import os
os.chdir(r"C:\Users\Yong\Desktop\ZMQ_phosopho")
os.chdir(os.getcwd())
f = open("pFind.protein", 'r')
b = open("report.txt", 'w')
tab = f.readlines()
b.write("\t".join([
    "ID", "Sequence", "Protein", "Modification", "Modi_site", "Spectra Number",
    "Final score"
]))
b.write("\n")

for i in range(2, len(tab)):
    if len(tab[i].strip().split("\t")) > 8 and tab[i].strip().split(
            "\t")[8] == "IARS/":  # import protein name
        peptide_start_site = int(tab[i].strip().split("\t")[9][:tab[
            i].strip().split("\t")[9].find(",")])
        modi = tab[i].strip().split("\t")[6]
        modi_list = modi.split(";")
        site_infor = []
        for mod in modi_list:
            if mod[mod.find(",") + 1:mod.find("[")] == "Phospho":  # modi name
                aa = mod[mod.find("[") + 1:mod.find("]")]
                site = int(mod[:mod.find(",")])
                site_infor.append(
                    "IARS[" + str(peptide_start_site + site) + "]")
        SITE_IFOR = ";".join(site_infor)
        peptide_ID = tab[i].strip().split("\t")[0]
        peptide_seq = tab[i].strip().split("\t")[1]
        final_score = tab[i].strip().split("\t")[5]
        peptide_mod = tab[i].strip().split("\t")[6]
        peptide_spec_num = tab[i].strip().split("\t")[-1]
        protein_name = tab[i].strip().split("\t")[8]
        if SITE_IFOR:
            b.write("\t".join([
                peptide_ID, peptide_seq, protein_name, peptide_mod, SITE_IFOR,
                peptide_spec_num, final_score
            ]))
            b.write("\n")
b.close()
f.close()

phospho_list = open("report.txt", 'r').readlines()
b = open("report.txt", 'a')
b.write("\t".join(["site name", "peptide_num", "site_spectra", "Best score"]))
b.write("\t".join("\n"))
all_site = []
for n in range(2, len(phospho_list)):
    site = phospho_list[n].strip().split("\t")[-3]
    if site.find(";"):
        site = site.split(";")
    all_site += site
print(all_site)
all_site = list(set(all_site))
all_site.sort()
print(all_site)

for name in all_site:
    spectra_num_list = []
    final_score_list = []
    for n in range(2, len(phospho_list)):
        site = phospho_list[n].strip().split("\t")[-3]
        if site.find(";") and name in site.split(";"):
            spectra_num_list.append(
                int(phospho_list[n].strip().split("\t")[-2]))
            final_score_list.append(
                float(phospho_list[n].strip().split("\t")[-1]))
        elif site == name:
            spectra_num_list.append(
                int(phospho_list[n].strip().split("\t")[-2]))
            final_score_list.append(
                float(phospho_list[n].strip().split("\t")[-1]))
    print(spectra_num_list)
    print(final_score_list)
    peptide_num = str(len(spectra_num_list))
    site_spectra = str(sum(spectra_num_list))
    Best_score = str(min(final_score_list))
    b.write("\t".join([name, peptide_num, site_spectra, Best_score]))
    b.write("\n")
b.close()
