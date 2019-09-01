# coding = utf-8

import os

#os.chdir("F:\MS_DATA_STORAGE\20190819")
flPath = r"F:\MS_DATA_STORAGE\20190819\multiNCE\BSA_DSSO_INCLU_RANDOM_REVERSE_R1.ms2"


def findMZinmzDic(mz, charge, mzdic):
    if mzdic == {}:
        return False
    else:
        if (mz, charge) in mzdic:
            return True
        else:
            mzList = sorted([x[0] for x in list(mzdic.keys())])
            #print(mzList)
            deltaMass = mz * 10 / 1000000
            lowMZ, upMZ = mz - deltaMass, mz + deltaMass
            if lowMZ > mzList[-1]:
                return False
            else:
                i = 0
                while i < len(mzList):
                    if upMZ < mzList[i]:
                        return False
                    elif lowMZ <= mzList[i]:
                        if (mzList[i], charge) in mzdic:
                            return True, mzList[i]
                        else:
                            return False
                    else:
                        i += 1
                if i == len(mzList):
                    return False


with open(flPath) as file:
    f = file.readlines()
repDic = {}
i = 0
while i < len(f):
    if f[i][0] != "S":
        i += 1
    else:
        monoMZ = float(f[i].split("\t")[-1].strip())
        #print(monoMZ)
        charge = f[i+9].split("\t")[1].strip()
        #print(monoMZ, charge)
        if not findMZinmzDic(monoMZ, charge, repDic):
            repDic[(monoMZ, charge)] = 1
        else:
            if findMZinmzDic(monoMZ, charge, repDic) == True:
                repDic[(monoMZ, charge)] += 1
            else:
                theroMZ= findMZinmzDic(monoMZ, charge, repDic)[-1]
                repDic[(theroMZ, charge)] += 1
        i += 9

print(len(repDic))