import os
os.chdir(r"E:\Script\cngp_modeling\20190520\DSSaRgo\local2refine\sd1169_3u28_relaxed_0176.topConf_NO1")
wantChains = ["C"]
#input_file = "3u28_relaxed.pdb"




def writeInfo(chain, pdb, flToWT):
    for line in pdb:
        if line[:3] == "TER" and str(line[21]) == chain:
            flToWT.write(line)
            break
        elif line[:4] == "ATOM" and str(line[21]) == chain:
            flToWT.write(line)
        else:
            continue


def main():
    for fl in os.listdir("./"):
        if fl[-6] == "7":
            input_file = fl
            pdb = open(input_file, 'r').readlines()
            b = open(input_file[:-4] + "_"+ "".join(wantChains) + ".pdb", 'w')

            for chain in wantChains:
                writeInfo(chain, pdb, b)

            b.close()


if __name__ == "__main__":
    main()