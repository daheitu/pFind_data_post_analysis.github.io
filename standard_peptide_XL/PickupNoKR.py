# coding = utf-8

import os

os.chdir(r"I:\CY\3_2_3_2")
b = open("filter.txt", 'w')
wList = []
with open("ans.txt") as fl:
    f = fl.readlines()
for line in f:
    if "K" not in line and 'R' not in line and "M" not in line:
        wList.append(line)
bList = sorted(wList, key = lambda x:float(x.strip().split("\t")[-1]), reverse = True)
for line in bList:
    b.write(line)
b.close()