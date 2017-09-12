import os
os.chdir(r"D:\softwareData\plink2\search_task_2017.09.05.15.20.16_BSA_low_concentration\reports")
filename = os.listdir(os.getcwd())
for fl in filename:
    if fl[-22:-4] == "cross-linked_sites":
        f=open(fl)
        b = open("report.txt", 'w')
        site_table=f.readlines()
        ordr_pos=[]
        sit_list=[]
        Peptide=[]
        for i in range(2,len(site_table)):
            if len([x for x in site_table[i].strip().split(',') if x !=''])==4:
                ordr_pos.append(i)
        ordr_pos.append(len(site_table))
        raw_name_list=[]
        for m in range(len(ordr_pos)-1):
            for k in range(ordr_pos[m]+1,ordr_pos[m+1]):
                a=site_table[k].strip().split(',')[2].find(".")
                raw_name_list.append(site_table[k].strip().split(',')[2][:a])
                raw_name_list = list(set(raw_name_list))
        raw_name_list.sort()
        print(raw_name_list)
        col = ["linked_site", "total_spec", "best_svm_score", "peptide", "peptide_mass", "prote_type"]
        for i in range(len(raw_name_list)):
            col.append(str(raw_name_list[i])+"_spec")
            col.append(str(raw_name_list[i])+"_svm")
            col.append(str(raw_name_list[i])+"unique_pep_num")

        b.write('\t'.join(col))
        b.write('\n')
        m = 0
        link_site_total_dic={}
        for l in range(len(ordr_pos) - 1):
            linked_site = site_table[ordr_pos[l]].strip().split(',')[1]
            link_site_total_dic[linked_site]= [
                site_table[ordr_pos[l]].strip().split(',')[1],
                site_table[ordr_pos[l]].strip().split(',')[3],
                site_table[ordr_pos[l] + 1].strip().split(',')[9],
                site_table[ordr_pos[l] + 1].strip().split(',')[5],
                site_table[ordr_pos[l] + 1].strip().split(',')[4]]
            m = linked_site.find("(")
            n = linked_site.find(")")
            p = linked_site.find("-")
            x = linked_site.find("(", p)
            y = linked_site.find(")", p)
            if linked_site[:m] == linked_site[p + 1:x]:
                link_site_total_dic[linked_site].append("intra")
            else:
                link_site_total_dic[linked_site].append("inter")
            raw_sub_spectra_dic={}
            raw_sub_svm_dic = {}
            raw_sub_peptide_dic = {}
            for raw in raw_name_list:
                raw_sub_spectra_dic[raw] = 0
                raw_sub_svm_dic[raw] = []
                raw_sub_peptide_dic[raw] = []
                for j in range(ordr_pos[l] + 1, ordr_pos[l + 1]):
                    if raw == site_table[j].strip().split(',')[2][:site_table[j].strip().split(',')[2].find(".")]:
                        raw_sub_spectra_dic[raw] += 1
                        raw_sub_svm_dic[raw].append(site_table[j].strip().split(',')[9])
                        raw_sub_peptide_dic[raw].append(site_table[j].strip().split(',')[5])
                if raw_sub_spectra_dic[raw] == 0:
                    link_site_total_dic[linked_site].append("")
                    link_site_total_dic[linked_site].append("")
                    link_site_total_dic[linked_site].append("")
                else:
                    link_site_total_dic[linked_site].append(str(raw_sub_spectra_dic[raw]))
                    link_site_total_dic[linked_site].append(raw_sub_svm_dic[raw][0])
                    link_site_total_dic[linked_site].append(str(len(list(set(raw_sub_peptide_dic[raw])))))
            #print(link_site_total_dic[linked_site])
            b.write('\t'.join(link_site_total_dic[linked_site]))
            b.write('\n')
        b.close()

        c = open("report.txt", 'a')
        rep_table = open("report.txt").readlines()
        print(rep_table)
        col_dic = {}
        total_spectra = 0
        total_colom = len(rep_table[0].strip().split("\t"))
        intra_num = 0
        for i in range(1, len(rep_table)):
            total_spectra += int(rep_table[i].strip("\n").split("\t")[1])
            if rep_table[i].strip("\n").split("\t")[5] == "intra":
                intra_num += 1

        col_dic[5] = float(intra_num) / (len(rep_table) - 1)
        col_dic[0] = len(rep_table) - 1
        col_dic[1] = total_spectra
        col_dic[2] = ""
        col_dic[3] = ""
        col_dic[4] = ""

        column_sub_dic={}
        for k in [6, 8]:
            for j in range(k, total_colom, 3):
                column_sub_dic[j]=[]
                for i in range(1, len(rep_table)):
                    if rep_table[i].strip("\n").split("\t")[j]:
                        column_sub_dic[j].append(int(rep_table[i].strip("\n").split("\t")[j]))
                col_dic[j] = sum(column_sub_dic[j])

        for j in range(7, total_colom, 3):
            column_sub_dic[j] = []
            for i in range(1, len(rep_table)):
                if rep_table[i].strip("\n").split("\t")[j]:
                    column_sub_dic[j].append(
                        round(float(rep_table[i].strip("\n").split("\t")[j]), 2))
            col_dic[j] = str(str(min(column_sub_dic[j])) + ',' + str(max(column_sub_dic[j])))
        last = []
        for k in range(total_colom):
            last.append(str(col_dic[k]))
        c.write("\t".join(last))
        print(last)
        c.close()