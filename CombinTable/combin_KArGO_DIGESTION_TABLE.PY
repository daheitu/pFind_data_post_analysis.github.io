# coding = utf8

import os

os.chdir(r'F:\OneDrive\kargo\TABLES\KArGO_DIGESTION')

f = open("KArGO_T.csv", 'r').readlines()
b = open("comTable.csv", 'w')

def extracInfo(flName):
    f = open(flName, 'r').readlines()
    infoDic = {}
    
    for line in f[1:]:
        if line.strip():
            lineList = line.strip().split(",")
            site = lineList[0]
            infoDic[site] = lineList
    return infoDic


KTdic = extracInfo("KArGO_T.csv")
KATdic = extracInfo("KArGO_AT.csv")
#print(KATdic)

combDic = {}
totalSite = sorted(list(set(KTdic.keys()) | set(KATdic.keys())))

print(len(totalSite))
for site in totalSite:
    if site in KTdic and site in KATdic:
        tSpec = int(KTdic[site][1]) + int(KATdic[site][1])
        mEvalue = min(float(KTdic[site][2]), float(KATdic[site][2]))
        mSvm = min(float(KTdic[site][3]), float(KATdic[site][3]))
        combDic[site] = [site, tSpec, mEvalue, mSvm]
        combDic[site].extend(KTdic[site][4:])
    elif site in KTdic:
        combDic[site] = KTdic[site]
    else:
        combDic[site] = KATdic[site]

for site in combDic:
    b.write(",".join([str(ele) for ele in combDic[site]])+"\n")

b.close()
