from urllib import request
import os
# from read_csv import get_target

# target = [45, 39, 47, 32, 35, 38, 37, 33, 40, 44]
ftpIP = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014337"
response = request.urlopen(ftpIP)
a = response.read().decode('utf-8').split("\n")
# print(a)
for line in a:
    name = line.split(" ")[-1].strip()
    if name.endswith(".raw"):# and "DSS_FAIMS_5060" in name:
        # if int(name[:-4].split('_')[-1]) in target:
        # if "re" not in name[-6:]:
        print(name)
        path = os.path.join(ftpIP, name)
        print(path)
        cmd = "wget " + path
        print(cmd)
        # os.system(cmd)
# ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD013928/JinLiang10366422_XL_SCX_v3_Fr_18.msf