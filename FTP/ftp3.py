import imp
from urllib import request
import os
import re
# from read_csv import get_target

# target = [45, 39, 47, 32, 35, 38, 37, 33, 40, 44]
ftpIP = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD014877/"
file_list = os.listdir("./")
response = request.urlopen(ftpIP)
a = response.read().decode('utf-8').split("\n")
# print(a)
for line in a:
    # name = line.split(" ")[-1].strip()
    if "href=" in line:
        print(line)
    # if name.endswith(".raw"):# and "Trypsin" in name and "GluC" not in file_list:
    #     # if int(name[:-4].split('_')[-1]) in target:
    #     # if "re" not in name[-6:]:
    #     print(name)
    #     path = os.path.join(ftpIP, name)
    #     print(path)
    #     cmd = "wget " + path
    #     print(cmd)
        # os.system(cmd)
# ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD013928/JinLiang10366422_XL_SCX_v3_Fr_18.msf
# https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD014877/20180920_QX3_JoMu_SA_LC12-7_uPAC200cm_HeLa+iRT_F3.raw