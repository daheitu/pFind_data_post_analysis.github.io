
import os
os.chdir(r"G:\20180724\AR_ARGO2")
filenames = os.listdir(os.getcwd())
comb = open("BSA_tryp_8_time.mgf", 'w')
for name in filenames:
    if name[-4:] == ".mgf":
        for line in open(name):
            comb.writelines(line)
comb.close
