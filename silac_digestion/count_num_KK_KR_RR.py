import os

os.chdir(r"E:\Script\xwalk_output\count")
linkerCutoffDic = {"KK_60" : 24, "KR_60" : 32.4, "RR_60" : 33.7}
disCufOff = 24


def countnum(flpath, disCufOff):
    f = open(flpath, 'r').readlines()
    if len(f[:-1]) == 0:
        return 0
    
    else:
        #print(f[0].split("\t"))
        tlnum = 0
        for line in f[:-1]:
            lineList = line.split("\t")
            disED = lineList[-5]
            disSADA = lineList[-4]
            if float(disED) < disCufOff:
                tlnum += 1
    return tlnum 


repDic = {}
for (root, dirs, files) in os.walk("./"):
    for fl in files:
        pdbName = fl[:4]
        linkType = root.split("\\")[-1][2:]
        print(linkType in linkerCutoffDic)
        if linkType in linkerCutoffDic:
            print(os.path.join(root, fl))
            tlnum = countnum(os.path.join(root, fl), linkerCutoffDic[linkType])
            if pdbName not in repDic:
                repDic[pdbName] = {}
                repDic[pdbName][linkType] = tlnum
            else:
                repDic[pdbName][linkType] = tlnum
print(repDic)
b = open("report_KK_KR_RR.txt", 'w')
b.write("\t".join(["pdb", "numKK", "numRR", "numKR"])+"\n")
for pdb in repDic:
    if "KK_60" in repDic[pdb]:
        numKK = repDic[pdb]["KK_60"]
    else:
        numKK = "No"
        #print(pdb)
    if "RR_60" in repDic[pdb]:
        numRR = repDic[pdb]["RR_60"]
    else:
        numRR = "No"
        print(pdb)
    if "KR_60"  in  repDic[pdb]:
        numKR = repDic[pdb]["KR_60"]
    else:
        numKR = "No"
        print(pdb)
    
    wlist = [pdb, numKK, numRR, numKR]
    
    b.write("\t".join([str(ele) for ele in wlist]) + "\n")

b.close()