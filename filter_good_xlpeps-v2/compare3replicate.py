# coding = utf-8

import os

path_r2 = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.26.14.53.44_CNGP_DSSO_preID_R2"
path_r3 = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.26.15.05.42_CNGP_DSSO_preID_R3"


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


def mergeInclison():
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


def extracOrgInfo(path):
    setIONS = set()
    flpath = os.path.join(path, "inclusion_list.csv")
    f = open(flpath, 'r').readlines()
    for line in f[1:]:
        lineList = line.split(",")
        pep = lineList[0]
        charge = lineList[1]
        mod = lineList[2]
        setIONS.add((pep, charge, mod))
    return setIONS, f


def mergeOrginclu():
    setR2, f2 = extracOrgInfo(path_r2)
    setR3, f3 = extracOrgInfo(path_r3)
    ov = setR2 & setR3
    set2Only = setR2 - ov
    b = open(os.path.join(path_r3, "merge_inclusion.csv"), 'w')
    for line in f3:
        b.write(line)
    
    for line in f2[1:]:
        lineList = line.split(",")
        pep = lineList[0]
        charge = lineList[1]
        mod = lineList[2]
        if (pep, charge, mod) in set2Only:
            b.write(line)
    b.close()


if __name__ == "__main__":
    mergeOrginclu()

