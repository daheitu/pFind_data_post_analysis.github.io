
import os

ms3_path =r"D:\wechat_files\WeChat Files\yunzaitianya2010\FileStorage\File\2022-09\peptide1N_1C_HCD.ms3" 
f = open(ms3_path).readlines()

b = open(os.path.join(os.path.dirname(ms3_path), os.path.basename(ms3_path)[:-4]+"ms3.mgf"), 'w')

i = 0
while i < len(f):
    if not f[i].startswith("S"):
        i += 1
    else:
        scan = str(int(f[i].split("\t")[1]))
        pre_mz = f[i].split("\t")[-1].strip()
        rt = f[i+1].split("\t")[-1].strip()
        charge = f[i+9].split("\t")[1]
        b.write("BEGIN IONS\n")
        b.write("TITLE=%s.%s.%s.%s.0.dta\n" % (os.path.basename(ms3_path)[:-4], scan, scan, charge))
        b.write("CHARGE=%s+\n" % charge)
        b.write("RTINSECONDS=%s\n" % rt)
        b.write("PEPMASS=%s\n" % pre_mz)
        p = i + 10
        while p < len(f):
            if f[p].startswith("S"):
                b.write("END IONS\n")
                break
            elif not f[p][0].isdigit():
                p += 1
            else:
                
                b.write(" ".join(f[p].split(" ")[:2])+"\n")
                p += 1
            
        i = p
    
b.write("END IONS\n")
b.close()