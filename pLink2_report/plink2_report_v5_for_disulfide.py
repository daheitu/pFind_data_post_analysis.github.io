# coding = utf-8

## created by Yong Cao


import os

os.chdir(r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.24.21.18.43_Integrin_LTG\reports")


def reportXLfile(openedXLfile, flTOwrt):
    f = openedXLfile; b = flTOwrt
    i = 2
    while i < len(f):
        lineList = f[i].split(",")
        if len(lineList) == 4:
            linkPair = [lineList[1]]
            if "SameSet" in f[i+1] or "SubSet" in f[i+1]:
                m = i + 1
                while m < len(f):
                    
                linkPair.append(f[i+1].split(",")[1])
                i += 1
            else:
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
                i = p
                site = "|".join(linkPair)
                wList = [site, numSpec, bestEvale, bestSVM, pep, spec1, spec2]
                b.write(",".join(wList) + "\n")
                
        else:
            print("wrong")


def main():
    for fl in os.listdir("./"):
        if fl.endswith("filtered_cross-linked_sites.csv"):
            xlfl = open(fl, 'r').readlines()
        if fl.endswith("filtered_loop-linked_sites.csv"):
            lpfl = open(fl, 'r').readlines()
    b = open("report.csv", 'w')
    b.write("\t".join(["linkPair", "numSpec", "bestEvale", "bestSVM", "pep", "spec1", "spec2"])+"\n")
    reportXLfile(xlfl, b)
    reportXLfile(lpfl, b)
    b.close()


if __name__ == "__main__":
    main()




