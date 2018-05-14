
import os
os.chdir(r"F:\A70_PEPTIDE\DK")
filenames = os.listdir(os.getcwd())
comb = open("combine.mgf", 'w')
for name in filenames:
    if name[-4:] == ".mgf":
        for line in open(name):
            comb.writelines(line)
comb.close


