# coding: utf-8

# author Yong CAO

import os


wk_dir = r'Z:\download_from_paper\mann_2020_nature'
demo_pparse = r"E:\pFindStudio\pLink2.3.9_0415\bin\pParse2021110111537.para"
os.chdir(r"E:\pFindStudio\pLink2.3.9_0415\bin")

# for root, dirs, fls in os.walk(wk_dir):
#     print(root)

def write_pparse(fl_path):
    b = open(os.path.join(fl_path, "pParse.para"), 'w')
    f = open(demo_pparse).readlines()
    for line in f:
        if line.startswith("datapath"):
            b.write("datapath = " + fl_path + "\n")
        else:
            b.write(line)
    b.close()



for fl in os.listdir(wk_dir):
    fl_path = os.path.join(wk_dir, fl)
    if os.path.isdir(fl_path):
        write_pparse(fl_path)
        parse_cmd = "pParse.exe " + os.path.join(fl_path, "pParse.para")
        os.system(parse_cmd)
