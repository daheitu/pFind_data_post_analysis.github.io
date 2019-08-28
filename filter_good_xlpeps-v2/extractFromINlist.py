# coding = utf-8

import os

os.chdir(r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.08.19.19.57.19")
specCF = 2
svmCF = 0.001
evalueCF = 0.001

b = open("inclusion.csv", 'w')
b.write("Compound,Formula,Adduct,m/z,z,t start (min),t stop (min)\n")


with open("./inclusion_list.csv", 'r') as file:
    f = file.readlines()
for line in f[1:]:
    lineList = line.strip().split(",")
    specNum = int(lineList[6])
    evalueBest = float(lineList[7])
    svmBest = float(lineList[8])
    pq_rt_min = float(lineList[10])
    pq_rt_max = float(lineList[11])
    pep = lineList[0]
    chg = lineList[1]
    mod = lineList[2]
    theo_mz = lineList[3]
    eluteTimeWin = pq_rt_max - pq_rt_min
    if specNum >= specCF and svmBest < svmCF and evalueBest < evalueCF:
        if eluteTimeWin > 40:
            eluateM = (pq_rt_max + pq_rt_min) / 2
            wList = [pep, "", "+H", theo_mz, chg, (eluateM - 300)/60, (eluateM+300)/60]
        else:
            eluateM = (pq_rt_max + pq_rt_min) / 2
            wList = [pep, "", "+H", theo_mz, chg, (eluateM - 300)/60, (eluateM+300)/60]
        #print(lineList[10], lineList[11])
        
        b.write(",".join([str(ele) for ele in wList]) + "\n")
b.close() 