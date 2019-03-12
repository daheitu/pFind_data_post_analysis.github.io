import os

os.chdir(r"E:\Script\test_treport ions\DSSO")

f = open('cy_1.spectra.xls', 'r').readlines()

print(f[2].strip().split("\t"))
scanTargetList = []
for line in f:
    if line[0].isdigit():
        lineList = line.strip().split("\t")
        title = lineList[1]
        scan = title.split('.')[1]
        scanTargetList.append(scan)
    else:
        continue

print(len(scanTargetList))

def writeTargetSpec(scanTargetList, orginal_mgf, flToWrite):
    mgf = open(orginal_mgf, 'r').readlines()
    i = 0
    while i < len(mgf):
        if mgf[i].strip() == "BEGIN IONS":
            begin_idx = i
        else:
            print("wrong")
        mix_num = mgf[i+1].split('.')[4]
        scan = mgf[i+1].split('.')[1]
        if scan in scanTargetList and mix_num == '0':
            p = i
            while p < len(mgf):
                if mgf[p].strip() == "END IONS":
                    flToWrite.write(mgf[p])
                    break
                else:
                    flToWrite.write(mgf[p])
                p += 1
        else:
            p = i + 5
            while p < len(mgf):
                if mgf[p].strip() == "END IONS":
                    break
                
                p += 1
        
        i = p + 1

rep = open('rep.mgf', 'w')
writeTargetSpec(scanTargetList, 'BSA_DSSO_BIO_R3_T2_HCDFT.mgf', rep)
rep.close()