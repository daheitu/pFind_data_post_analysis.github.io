import os
os.chdir(r"D:\softwareData\plink2\20160703\search_task_2016.08.03.09.16.26_BSA_BH\reports")
f=open("bsa_2016.08.03.filtered_cross-linked_spectra.csv")
inter_table=f.readlines()
inter_scan_numb=[]
b=open(r"F:\20160719\BSA\hepes_borate\BSA_argo2_hepes(BORATE)_50mm.ms1")
ms1=b.readlines()
for i in range(1,len(inter_table)):
    if len(inter_table[i].strip().split(","))>3 and inter_table[i].strip().split(",")[-3]=="1":
        scan_num=int(inter_table[i].strip().split(",")[-2])
        charge=inter_table[i].strip().split(",")[2]
        moverz=(float(inter_table[i].strip().split(",")[3])+int(charge))/int(charge)
        for k in range(200000,len(ms1)):
            if ms1[k][0]=="S" and int(ms1[k].strip().split(" ")[-1])>scan_num-15 and int(ms1[k].strip().split(" ")[-1])<scan_num :

            else:
                if
print inter_scan_numb.sort()
print inter_scan_numb


