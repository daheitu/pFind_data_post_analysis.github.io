import os

os.chdir(r"D:\softwareData\plink2\pLink_task_2018.03.05.07.48.38_KARGO_T\reports")

f = open("ALPK1_con_2018.03.05.filtered_mono-linked_sites.csv")


def site_list_process(site_list):
    i = 0
    while i < len(site_list):
        if "REVERSE" in site_list[i] or "gi|CON" in site_list[i]:
            site_list.remove(site_list[i])
        else:
            i += 1
    return site_list


site_tgt_list = ["ALPK1(38)", "ALPK1(132)", "ALPK1(329)", "ALPK1(383)"]
n = 2
while n < len(f):
    line = f[n]
    line_list = line.rstrip("\n").split(",")
    if line_list[0].isdigit():
        sites = [line_list[1]]
        total_spec = line_list[3]         
        p = n + 1
        if f[p].rstrip("\n").split(",")[0] in ["SameSet", "SubSet"]:         
            while p < len(f) and f[p].rstrip("\n").split(",")[0] in [
                    "SameSet", "SubSet"]:
                if f[p].rstrip("\n").split(",")[0] == "SameSet":
                    sites.append(f[p].rstrip("\n").split(",")[1])
                else:
                    pass
                    p += 1
        else:
            pass
        sites = site_list_process(sites)
        in_tgt = False
        for site in sites:
            if site in site_tgt_list:
                in_tgt = True
                break                
            else:
                continue
        if in_tgt:
            while p < len(f) and f[p].rstrip("\n").split(",")[0] == "":
                

        else:
            while p < len(f) and not f[p].rstrip("\n").split(",")[0].isdigit():
                p += 1
        
    else:
        print(n)
    
    if f[p].rstrip("\n").split(",")[0] in ["SameSet", "SubSet"]:
        while p < len(f) and f[p].rstrip("\n").split(",")[0] in [
                "SameSet", "SubSet"
        ]:
            if f[p].rstrip("\n").split(",")[0] == "SameSet":
                site_list.append(f[p].rstrip("\n").split(",")[1])
            else:
                pass
            p += 1
    else:
        pass
    site = site_list_process(site_list)[0]
    if site == "":
        while p < len(f) and f[p].rstrip("\n").split(",")[0] == "":
            p += 1
    else:
        link_type = site_list_process(site_list)[1]
        best_svm_score = f[p].rstrip("\n").split(",")[9]
        pep_std = f[p].rstrip("\n").split(",")[5]
        while p < len(f) and f[p].rstrip("\n").split(",")[0] == "":
            line_list = f[p].rstrip("\n").split(",")
            raw_name = line_list[2][:line_list[2].find(".")]
            pep = line_list[5]
            evalue = line_list[8]

            if raw_name not in spec_dic:
                spec_dic[raw_name] = 1
            else:
                spec_dic[raw_name] += 1

            if raw_name not in pep_dic:
                pep_dic[raw_name] = [pep]
            else:
                if pep not in pep_dic[raw_name]:
                    pep_dic[raw_name].append(pep)
                else:
                    pass

            if raw_name not in evalue_dic:
                evalue_dic[raw_name] = float(evalue)
            else:
                if float(evalue) < evalue_dic[raw_name]:
                    evalue_dic[raw_name] = float(evalue)
                else:
                    pass
            p += 1

        min_evalue = 1
        for key in evalue_dic:
            if evalue_dic[key] < min_evalue:
                min_evalue = evalue_dic[key]
            else:
                continue

        rep_dic = {}
        for raw_name in raw_name_list:
            if raw_name not in spec_dic:
                rep_dic[raw_name] = ["", "", ""]
            else:
                rep_dic[raw_name] = [
                    str(spec_dic[raw_name]),
                    str(evalue_dic[raw_name]),
                    str(len(pep_dic[raw_name]))
                ]
        rep_list = [
            site, total_spec,
            str(min_evalue), best_svm_score, pep_std, link_type
        ]
        for raw_name in raw_name_list:
            rep_list.extend(rep_dic[raw_name])

        if float(rep_list[1]) > spec_cutoff and float(
                rep_list[2]) < E_value_cutoff:
            b.write("\t".join(rep_list))
            b.write("\n")
        else:
            pass

    n = p
b.close()
