# coding = utf-8

import os

path_r2 = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.11.23.14.08.46_CNGP_DSSO_preID_R1"
path_r3 = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.11.23.14.16.25_CNGP_DSSO_preID_R2"


def extractIncluList(path):
    specCF = 1
    svmCF = 0.001
    evalueCF = 1e-6

    b = open(os.path.join(path, "inclusion.csv"), 'w')
    b.write("Compound,Formula,Adduct,m/z,z,t start (min),t stop (min)\n")

    with open(os.path.join(path, "inclusion_list.csv"), 'r') as file:
        f = file.readlines()
    for line in f[1:]:
        lineList = line.strip().split(",")
        specNum = int(lineList[6])
        evalueBest = float(lineList[7])
        svmBest = float(lineList[8])
        #pq_rt_min = float(lineList[10])
        #pq_rt_max = float(lineList[11])
        fitChamoMean = float(lineList[12])
        pep = lineList[0]
        chg = lineList[1]
        #mod = lineList[2]
        theo_mz = lineList[3]
        if specNum >= specCF:
            if svmBest < svmCF or evalueBest < evalueCF:
                eluateM = fitChamoMean  # (pq_rt_max + pq_rt_min) / 2
                wList = [pep, "", "+H", theo_mz, chg, (eluateM - 120)/60, (eluateM+120)/60]
                b.write(",".join([str(ele) for ele in wList]) + "\n")
    b.close() 


def extractInfo(path):
    extractIncluList(path)
    pepList = []
    flpath = os.path.join(path, "inclusion.csv")
    f = open(flpath).readlines()
    for line in f[1:]:
        lineList = line.split(",")
        pep = lineList[0]
        thero_mz = lineList[3]
        pepList.append((pep, thero_mz))

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
    mergeInclison()
    mergeOrginclu()