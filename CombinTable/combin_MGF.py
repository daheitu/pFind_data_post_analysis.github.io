# coding = utf-8
import os
os.chdir(r"G:\20190402") # important; the path of .mgf file 
def combineMGF(flToW):
    flList = os.listdir(os.getcwd())
    for fl in flList:
        if fl[-8:] == "CDFT.mgf":
            print("The current mgf file is " + fl)
            for line in open(fl):
                flToW.writelines(line)
    


def main():
    combFLname = "bsa.mgf" # important name of combined file
    b = open(combFLname, 'w')
    combineMGF(b)
    b.close()


if __name__ == "__main__":
    main()