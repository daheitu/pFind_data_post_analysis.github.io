# coding = utf-8

## created by Yong Cao


import os

os.chdir(r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.07.01.15.13.27\reports")


def reportXLfile(openedXLfile, flTOwrt):
    f = openedXLfile; b = flTOwrt
    i = 2
    while i < len(f):
        lineList = f[i].split(",")
        if len(lineList) == 4:
            linkPair = lineList[1]
            numSpec = lineList[-1].strip()
            bestSVM = f[i+1].split(",")[9]
            pep = f[i+1].split(",")[5]
            spec1 = f[i+1].split(",")[2]
            if int(numSpec) > 1:
                spec2 = f[i+2].split(",")[2]
            else:
                spec2 = ""
            evalueList = []
            p = i + 1
            while p < len(f):
                if len(f[p].split(",")) == 4:
                    break
                else:
                    evalue = float(f[p].split(",")[8])
                    evalueList.append(evalue)
                    p += 1
            bestEvale = str(min(evalueList))
            wList = [linkPair, numSpec, bestEvale, bestSVM, pep, spec1, spec2]
            b.write("\t".join(wList) + "\n")
            i = p
        else:
            print("wrong")


def main():
    xlfl = open(r"Integrin_2019.07.01.filtered_cross-linked_sites.csv", 'r').readlines()
    lpfl = open(r"Integrin_2019.07.01.filtered_loop-linked_sites.csv", 'r').readlines()
    b = open("report.txt", 'w')
    b.write("\t".join(["linkPair", "numSpec", "bestEvale", "bestSVM", "pep", "spec1", "spec2"])+"\n")
    reportXLfile(xlfl, b)
    reportXLfile(lpfl, b)
    b.close()


if __name__ == "__main__":
    main()




