import os

scoreCutoff = 0.5
os.chdir(r"C:\Users\Yong\Documents\pLink\pLink_task_2019.03.13.11.51.06\reports")
mgfPath = "G:\DSSO_1008\hechengtiaduan"

def getRetantionTime(IDSpec, mgfPath):
    raw_scan_dic = {}
    raw_rt_dic = {}
    for line in IDSpec[1:]:
        line_list = line.strip().split(",")
        title = line_list[1]
        e_value = float(line_list[9])
        score = float(line_list[10])
        raw_name, scan_num = title.split(".")[:2]
        raw_name += "_HCDFT.mgf"
        if  score < scoreCutoff:
            if raw_name not in raw_scan_dic:
                raw_scan_dic[raw_name] = [int(scan_num)]
            else:
                if scan_num not in raw_scan_dic[raw_name]:
                    raw_scan_dic[raw_name].append(int(scan_num))
                else:
                    pass
        else:
            continue


    raw_list = sorted(list(raw_scan_dic.keys()))

    for raw in raw_list:
        print(raw)
        rep = open(raw + "_rt.txt", "w")
        rep.write(raw + "\n")
        raw_scan_dic[raw].sort()
        raw_path = os.path.join(mgfPath, raw)
        raw_rt_dic[raw] = []
        mgf = open(raw_path, "r").readlines()
        scan_list = sorted(raw_scan_dic[raw])
        m = 0; n = 0
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
                    rep.write(str(trg_pos)+"\t"+str(rt_Sec) + "\n" )
            else:
                m += 1
        
        if len(raw_scan_dic[raw]) != len(raw_rt_dic[raw]):
            print("this raw have problem %s" % raw)
        else:
            continue    


for fl in os.listdir(r"./"):
    if fl[-24:] == "cross-linked_spectra.csv":
        f = open(fl, 'r').readlines()
        getRetantionTime(f, mgfPath)
    else:
        continue
                    