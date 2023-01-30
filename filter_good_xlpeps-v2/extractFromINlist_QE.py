# coding = utf-8

import os

os.chdir(r"K:\20200820_PRM\BSA\bsa_output_plink_dsso")
#C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.26.14.53.44_CNGP_DSSO_preID_R2
#C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.26.15.05.42_CNGP_DSSO_preID_R3
specCF = 1
svmCF = 1
evalueCF = 1

# nceList = list(range(20, 40, 2))

b = open("inclusion.csv", 'w')
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
    theo_mz = round(float(lineList[3]), 4)
    eluteTimeWin = pq_rt_max - pq_rt_min
    if specNum >= specCF and svmBest < svmCF and evalueBest <= evalueCF:
        inc_list = [""] * 12
        inc_list[0] = str(theo_mz)
        inc_list[4] = chg
        inc_list[5] = "Positive"
        inc_list[6] = str((fitChamoMean - 150)/60)
        inc_list[7] = str((fitChamoMean+150)/60)
        inc_list[8] = "" # 能量
        inc_list[9] = ""

        # wList = [str(theo_mz), "", "",  "", chg, "Positive", (fitChamoMean - 150)/60, (fitChamoMean+150)/60, "", "NCE",""]
        # for nce in nceList:
        #     wList[-3] = nce
        b.write(",".join([str(ele) for ele in inc_list]) + "\n")
     
b.close() 
