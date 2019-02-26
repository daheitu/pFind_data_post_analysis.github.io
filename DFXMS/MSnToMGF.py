import os

os.chdir(r"G:\DSSO_1008\CID_BASED_Methods\test")

f = open("DSSO_CID_MS2_CID_MS3_T1.ms2", 'r').readlines()
raw_name = 'DSSO_CID_MS2_CID_MS3_T1'
b = open(raw_name + ".mgf", 'w')
print(f[16].strip().split("\t"))
i= 16
while i < len(f):
    print(i)
    if f[i][0] == "S":
        scan = f[i].strip().split("\t")[1]
        precu = f[i].strip().split("\t")[-1]
    else:
        print('wrong')
    
    chrg_list = []
    m = i
    while m < len(f):
        if f[m][0].isdigit():
            break
        else:
            if f[m][:9] == "I	RetTime":
                ret_time = str(float(f[m].strip().split("\t")[-1]) * 60)
            elif f[m][0] == "Z":
                chrg = f[m].split("\t")[1]
                chrg_list.append(chrg)
        
        m += 1
    
    p = m
    if len(chrg_list) == 1:
        chrg = chrg_list[0]
    else:
        chrg = chrg_list[0]
        print(scan, len(chrg_list))
    title = ".".join([raw_name, scan, scan, chrg, "0.dta"])
    b.write("BEGIN IONS" + "\n")
    b.write("TITLE=" + title + "\n")
    b.write("CHARGE=" + chrg + "+\n")
    b.write("RTINSECONDS=" +  ret_time + '\n')
    b.write("PEPMASS="  + precu + '\n')
    while p < len(f):
        if f[p][0] == "S":
            break
        else:
            b.write(f[p])
        p += 1
    
    i = p
    b.write("END IONS\n")
b.close()


"""
for line in f[16:25]:
    print(line[:-1])
    """