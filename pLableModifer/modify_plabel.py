# coding = utf-8

import os, shutil

os.chdir(r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.05.07.10.27.20MG_LTG\tmps\complex_ss0")

def reNamePath(plabelFL, realpath):
    f = open(plabelFL, 'r').readlines()
    repName = plabelFL[:-7]+"re.plabel"
    b = open(repName, 'w')
    for line in f:
        if line[:9] == "File_Path":
            mgfName = line.split("\\")[-1]
            b.write("File_Path=" + realpath + "\\" + mgfName)
        else:
            b.write(line)


def moveplabel():
    os.mkdir("plabel")
    for fl in os.listdir("./"):
        if fl[-6:] == "plabel" and fl[-9:-7] != "re":
            shutil.move(fl, "./plabel")
        else:
            pass


def main():
    flList = os.listdir("./")
    realPath = r"E:\HGMG_DATA"
    for fl in flList:
        if fl[-6:] == "plabel":
            reNamePath(fl, realPath)
        else:
            pass
    moveplabel()


if __name__ == "__main__":
    main()