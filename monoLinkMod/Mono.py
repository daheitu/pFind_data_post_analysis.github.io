import os

os.chdir(r"D:\softwareData\plink2\pLink_task_2018.03.05.07.48.38_KARGO_T\reports")

f = open("ALPK1_con_2018.03.05.filtered_mono-linked_sites.csv", 'r').readlines()


def site_list_process(site_list):
    i = 0
    while i < len(site_list):
        if "REVERSE" in site_list[i] or "gi|CON" in site_list[i]:
            site_list.remove(site_list[i])
        else:
            i += 1
    return site_list


site_tgt_list = ["ALPK1(38)", "ALPK1(132)", "ALPK1(329)", "ALPK1(383)", "ALPK1(348)", "ALPK1(40)"]
raw_list = ["NC_KARGO", "N_KARGO"]
final_dic = {}
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
                w_site = site                
                print(w_site)
                in_tgt = True
                break                
            else:
                continue
        
        if in_tgt:
            print(p)
            final_dic[w_site] = {}
            while p < len(f) and f[p].rstrip("\n").split(",")[0] == "":
                in_raw = False
                for raw in raw_list:
                    if raw in f[p]:
                        in_raw = True
                        if raw not in final_dic[w_site]:
                            final_dic[w_site][raw] = 1
                        else:
                            final_dic[w_site][raw] += 1
                    else:
                        continue
                if not in_raw:
                    print(p)
                p += 1
            print(p)

        else:
            while p < len(f) and f[p].rstrip("\n").split(",")[0] == "":
                p += 1
        
    else:
        print(n)
    
    n = p

print(final_dic)
