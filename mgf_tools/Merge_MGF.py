# coding = utf-8
import os
os.chdir(r"M:\synthetic_pepteide_rawdata_NCE30\cov_cat\sensitivity_test\70_60") # important; the path of .mgf file 


def combineMGF(flToW):
    flList = os.listdir(os.getcwd())
    for fl in flList:
        if fl.endswith(".mgf"):
            print("The current mgf file is " + fl)
            for line in open(fl):
                flToW.writelines(line)


def main():
    combFLname = "70_60.mgf" # important name of combined file
    b = open(combFLname, 'w')
    combineMGF(b)
    b.close()


if __name__ == "__main__":
    main()