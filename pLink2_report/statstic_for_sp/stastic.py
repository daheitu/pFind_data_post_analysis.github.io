# coding = utf8

import os

os.chdir(r"/Users/yong/PycharmProjects/reports")

f = open("AR7_con_2019.07.03_Trypsin_DSSv4.txt"， 'r').readlines()
for line in f:
    print