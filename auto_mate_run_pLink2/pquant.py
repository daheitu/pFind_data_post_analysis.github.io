import os

plink_bin_path = r"E:\pFindStudio\pLink2.3.9_0415\bin"

os.chdir(plink_bin_path)
f = open('./pQuant_cfg.txt').readlines()
b = open('pquant_para.txt', 'w')
for line in f:
    if line.startswith("PATH_INI_ELEMENT"):
        b.write()