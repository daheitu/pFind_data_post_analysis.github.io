# coding = utf-8
# author: Yong CAO

import os
from filter_modification_comman_used import main

wkdir = r"E:\workspace"
protein = "sp|P01112|RASH_HUMAN"
modification_list = ["F[Y]", "C15H23F[C]", "Farnesyl[C]"]

for dr in os.listdir(wkdir):
    if dr.startswith("pFindTask111") and "open_modi" in dr:
        print(dr)
        t_path = os.path.join(wkdir, dr, "result")
        for mod in modification_list:
            main(t_path, protein, mod)
