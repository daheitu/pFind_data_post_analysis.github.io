# coding = utf-8

import os

os.chdir(r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.26.14.53.44_CNGP_DSSO_preID_R2")
#C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.26.14.53.44_CNGP_DSSO_preID_R2
#C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.26.15.05.42_CNGP_DSSO_preID_R3
specCF = 2
svmCF = 0.005
evalueCF = 0.01

nceList = list(range(20, 40, 2))

b = open("inclusion_QE.csv", 'w')
b.write("Mass [m/z],Formula [M],Formula type,Species,CS [z],Polarity,Start [min],End [min],(N)CE,(N)CE type,MSX ID,Comment\n")


with open("./inclusion_list.csv", 'r') as file:
    f = file.readlines()
for line in f[1:]:
    lineList = line.strip().split(",")
    specNum = int(lineList[6])
    evalueBest = float(lineList[7])
    svmBest = float(lineList[8])
    pq_rt_min = float(lineList[10])
    pq_rt_max = float(lineList[11])
    fitChamoMean = float(lineList[12])
    pep = lineList[0]
    chg = lineList[1]
    mod = lineList[2]
    theo_mz = lineList[3]
    eluteTimeWin = pq_rt_max - pq_rt_min
    if specNum >= specCF and svmBest < svmCF and evalueBest < evalueCF:
        wList = [theo_mz, "", "+H",  "", chg, "Positive", (fitChamoMean - 120)/60, (fitChamoMean+120)/60, "", "NCE",""]
        for nce in nceList:
            wList[-3] = nce
            b.write(",".join([str(ele) for ele in wList]) + "\n")
     
b.close() 
