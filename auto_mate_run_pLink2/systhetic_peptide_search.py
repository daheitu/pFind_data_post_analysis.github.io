import os

wd_dir = r"G:\msData\synthetic_pepteide_rawdata"


def remove_folder(path):
    # while os.path.exists(path):
    for root, dirs, fls in os.walk(path, topdown=False):
        for fl in fls:
            os.remove(os.path.join(root, fl))
        for dr in dirs:
            os.rmdir(os.path.join(root, dr))
    os.rmdir(path)


# remove_folder(r"G:\msData\synthetic_pepteide_rawdata\CV1_CF4_6\DSSO\CV1_CF4_6_DSSO")

# """
for root, dirs, fls in os.walk(wd_dir):
    linker = root.split("\\")[-1]
    if linker =="DSSO":
        for dr in dirs:
            if dr.endswith("DSSO"):
                dr_path = os.path.join(root, dr)
                print("We will remove the folder %s" % dr_path)
                remove_folder(dr_path)
        # print(root, dirs)
        # """