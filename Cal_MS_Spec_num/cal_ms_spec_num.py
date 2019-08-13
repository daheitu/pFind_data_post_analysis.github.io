# coding = utf8

import os

os.chdir(r"F:\MS_DATA_STORAGE\20190715_lumos\raw") # important ! change the path where are your .ms1, ms2

def main():
    flList = os.listdir("./")
    for fl in flList:
        if ".ms" in fl:
            rawName = fl[:-4]
            ms = fl[-3:]
            specCount = 0
            with open(fl) as file:
                f = file.readlines()
            for line in f:
                if line[0] == "S":
                    specCount += 1
            print("the Spec Num of " + rawName + " "+  ms + " is " + str(specCount))


if __name__ == "__main__":
    main()            