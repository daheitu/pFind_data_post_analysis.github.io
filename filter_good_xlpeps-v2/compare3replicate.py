# coding = utf-8

import os

path_r2 = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.10.02.18.29.48_Lacto_DSSO_R2"
path_r3 = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.10.02.18.35.44_Lacto_DSSO_R3"


def extractInfo(path):
    pepList = []
    flpath = os.path.join(path, "inclusion.csv")
    f = open(flpath).readlines()
    for line in f[1:]:
        lineList = line.split(",")
        pep = lineList[0]
        mz = lineList[3]
        pepList.append((pep, mz))

    return set(pepList)

setR2 = extractInfo(path_r2)
setR3 = extractInfo(path_r3)
ov = setR2 & setR3
set2Only = setR2 - ov
set3Only = setR3 - ov
print(len(ov), len(set2Only), len(set3Only))

b = open(os.path.join(path_r3, "inclusion.csv"), "a")

f = open(os.path.join(path_r2, "inclusion.csv")).readlines()
for line in f[1:]:
    lineList = line.split(",")
    pep = lineList[0]
    mz = lineList[3]
    if (pep, mz) in set2Only:
        b.write(line)

b.close()
