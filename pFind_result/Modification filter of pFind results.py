import os

wk_dir = r"E:\workspace\pFindTask97\result"
pro_looking = "hp_pro"  # important! enter the protein name you want search
modify = "hp3m[C]"  # important! enter the modification type you want search

rep_file_name = modify + "_" + pro_looking+"_filter.txt"


def filter_modification_from_pfind(path, tab, b):
    for line in tab[2:]:
        line_list = line.rstrip("\n").split("\t")
        if line_list[0].isdigit():
            if line_list[1] == pro_looking:
                print("The protein you are looking for in this list")
                pro_pos = tab.index(line)
                break
        else:
            continue
    print(pro_pos)
    for i in range(pro_pos + 1, len(tab)):
        line_list = tab[i].rstrip("\n").split("\t")
        if line_list[1] == "SubSet" or line_list[1] == "SameSet":
            continue
        elif line_list[0].isdigit():
            print(i)
            break
        else:
            protein_list = line_list[10].split("/")[:-1]
            if pro_looking not in protein_list:
                print("wrong")
            else:
                pro_index = protein_list.index(pro_looking)

            peptide_start_site_list = line_list[11].split("/")[:-1]
            peptide_start_site = int(
                peptide_start_site_list[pro_index].split(",")[0])

            site_infor = []
            modi_list = line_list[8].split(";")[:-1]
            for mod in modi_list:
                mod_type = mod.split(",")[1][:mod.split(",")[1].find("[")]
                if modify == mod_type:
                    print()  # modi name
                    aa = mod[mod.find("[") + 1:mod.find("]")]
                    site = int(mod.split(",")[0])
                    site_infor.append(pro_looking + "[" +
                                      str(peptide_start_site + site) + "]")
            if site_infor:
                SITE_IFOR = ";".join(site_infor)
                peptide_ID = line_list[2]
                peptide_seq = line_list[3]
                final_score = line_list[7]
                peptide_mod = line_list[8]
                peptide_spec_num = line_list[-1]
                protein_name = line_list[10]
                spectra_title = line_list[-3]
                b.write("\t".join([
                    peptide_ID, peptide_seq, protein_name, peptide_mod,
                    SITE_IFOR, peptide_spec_num, final_score, spectra_title
                ]))
                b.write("\n")
            else:
                continue
    return


def statistic_to_site_level():
    b = open("report.txt", 'r')
    phospho_list = b.readlines()
    b.close()
    b = open("report.txt", 'a')
    b.write("\t".join(
        ["site name", "peptide_num", "site_spectra", "Best score"]))
    b.write("\t".join("\n"))
    all_site = []
    for n in range(2, len(phospho_list)):
        site = phospho_list[n].strip().split("\t")[-4]
        if site.find(";"):
            site = site.split(";")
        all_site += site

    all_site = list(set(all_site))
    all_site.sort()

    for name in all_site:
        spectra_num_list = []
        final_score_list = []
        for n in range(2, len(phospho_list)):
            site = phospho_list[n].strip().split("\t")[-4]
            if site.find(";") and name in site.split(";"):
                spectra_num_list.append(
                    int(phospho_list[n].strip().split("\t")[-3]))
                final_score_list.append(
                    float(phospho_list[n].strip().split("\t")[-2]))
            elif site == name:
                spectra_num_list.append(
                    int(phospho_list[n].strip().split("\t")[-3]))
                final_score_list.append(
                    float(phospho_list[n].strip().split("\t")[-2]))
        print(spectra_num_list)
        print(final_score_list)
        peptide_num = str(len(spectra_num_list))
        site_spectra = str(sum(spectra_num_list))
        Best_score = str(min(final_score_list))
        b.write("\t".join([name, peptide_num, site_spectra, Best_score]))
        b.write("\n")
    b.close()
    return


def main():
    os.chdir(wk_dir)
    tab = open("pFind.protein", 'r').readlines()
    b = open(rep_file_name, 'w')
    b.write("\t".join([
        "ID", "Sequence", "Protein", "Modification", "Modi_site",
        "Spectra Number", "Final score", "Spectra title"
    ]))
    b.write("\n")
    path = os.getcwd()
    filter_modification_from_pfind(path, tab, b)
    statistic_to_site_level()
    b.close()


if __name__ == "__main__":
    main()
