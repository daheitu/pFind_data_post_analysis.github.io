import os
import shutil

aa_dir_dic = {"asp": "KD_60", "cys": "KC_60", "glu":"KE_60"}


os.chdir(r"E:\Script\xwalk_output\result")
targetd_path = r"E:\Script\xwalk_output"
for (root, dirpath, files) in os.walk("./"):
    for fl in files:
        for aa in aa_dir_dic:
            if aa in root:
                print(os.path.join(root, fl))
                shutil.move(os.path.join(root, fl), os.path.join(targetd_path, aa_dir_dic[aa]))
            