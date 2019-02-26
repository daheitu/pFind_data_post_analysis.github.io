import os


os.chdir(r"G:\DSSO_1008\DSS_PLINK1\2.report\sample1")

f = open("cy.spectra.xls", 'r').readlines()
raw_scan_dic = {}
raw_rt_dic = {}
i = 2
while i < len(f):
    line_list = f[i].strip().split("\t")
    title = line_list[1]
    e_value = line_list[5]
    dot_pos = title.find('.')
    dot_pos_next = title.find('.', dot_pos+1)
    raw_name = title[:dot_pos] + "_HCDFT.mgf"
    scan_num = title[dot_pos+1: dot_pos_next]
    if float(e_value) < 0.001:
        if raw_name not in raw_scan_dic:
            raw_scan_dic[raw_name] = [int(scan_num)]
        else:
            if scan_num not in raw_scan_dic[raw_name]:
                raw_scan_dic[raw_name].append(int(scan_num))
            else:
                pass
    else:
        pass
    
    i += 2

raw_list = list(raw_scan_dic.keys())

for raw in raw_list:
    # print(raw)
    rep = open(raw + "_rt.txt", "w")
    rep.write(raw + "\n")
    raw_scan_dic[raw].sort()
    raw_path = r'G:\DSSO_1008\%s'  % raw
    raw_rt_dic[raw] = []
    mgf = open(raw_path, "r").readlines()
    scan_list = sorted(raw_scan_dic[raw])
    m = 0
    n = 0
    while m < len(mgf) and n < len(scan_list):
        trg_pos = scan_list[n]
        if mgf[m].strip() == "BEGIN IONS" and mgf[m+1].strip()[-5:] == "0.dta":
            scan_num = int(mgf[m+1].split(".")[1])
            rt_Sec = round(float(mgf[m+3].strip().split("=")[1])/60, 2)
            if scan_num < trg_pos:
                m += 1
            elif scan_num > trg_pos:
                n += 1
            else:
                m += 1
                n += 1
                raw_rt_dic[raw].append(rt_Sec)
                rep.write(str(rt_Sec) + "\n" )
                # print(trg_pos)
        else:
            m += 1
    
    if len(raw_scan_dic[raw]) != len(raw_rt_dic[raw]):
        print(raw)
    else:
        continue    
            
       

    
    

