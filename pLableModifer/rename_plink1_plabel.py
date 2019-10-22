import os

os.chdir(r"G:\DSSO_1008\DSSO_PLINK1\2.report\sample1")


def rename_plabel(filename):
    name = filename[:filename.find(".plabel")]
    return name + "_re.plabel"


file_list = os.listdir(os.getcwd())
for file in file_list:
    if file[-13:] == "_inter.plabel":
        f = open(file, 'r').readlines()
        b = open(rename_plabel(file), 'w')
        for line in f:
            if line[:9] == "File_Path":
                # "File_Path=G:\DSSO_1008" + line[64:]
                b.write("File_Path=G:\DSSO_1008\\" + line[62:])
            else:
                b.write(line)

print("done")