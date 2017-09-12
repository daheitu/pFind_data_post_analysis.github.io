import os
os.chdir(r"E:\20170310\T3")
filenames = os.listdir(os.getcwd())
comb = open("T3.mgf",'w')
for name in filenames:
    if name[-4:]== ".mgf":
        for line in open(name):
            comb.writelines(line)
comb.close

