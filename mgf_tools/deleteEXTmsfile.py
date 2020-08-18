import os

path = r"F:\disulfide_Lumos"

for root, dirs, files in os.walk(path):
    for file in files:
        ends = file[file.find("."):]
        if ends in [".ms1", ".ms2", ".pf1", ".pf2", ".mgf", ".pf1idx", ".pf2idx", ".xtract"]:
            os.remove(os.path.join(root, file))