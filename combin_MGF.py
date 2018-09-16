
import os
os.chdir(r"F:\ALL ArGO data\ArGO\optimize condition\Buffer screen\Aldolase")
filenames = os.listdir(os.getcwd())
comb = open("BSA_buffer_Screen.mgf", 'w')
for name in filenames:
    if name[-4:] == ".mgf":
        for line in open(name):
            comb.writelines(line)
comb.close
