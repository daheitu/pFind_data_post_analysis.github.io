import os
import shutil

rep_path = r"G:\msData\SLu"

for root, dirs, fls in os.walk(rep_path, topdown= False):
    root_list = root.split("\\")
    if len(root_list) > 5 and root_list[4] == "pParse_para":
        # print(root)
        if root_list[-1] == "LuS_NP_output_BS3_heavy":
            print(root)
            shutil.rmtree(root)