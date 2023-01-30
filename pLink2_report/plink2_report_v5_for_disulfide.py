# coding = utf-8

## created by Yong Cao, 20190926


import os

os.chdir(r"C:\Users\Yong Cao\Documents\pLink\na_1013_TA\reports")


def continueLines(starPos, f, numSpec):
    bestSVM = f[starPos].split(",")[9]
    pep = f[starPos].split(",")[5]
    spec1 = f[starPos].split(",")[2]
    if int(numSpec) > 1:
        spec2 = f[starPos + 1].split(",")[2]
    else:
        spec2 = ""
    evalueList = []
    p = starPos
    while p < len(f):
        if len(f[p].split(",")) == 4:
            break
        else:
            evalue = float(f[p].split(",")[8])
            evalueList.append(evalue)
            p += 1
    bestEvale = str(min(evalueList))
    endPos = p
    return bestSVM, bestEvale, pep, spec1, spec2, endPos


def removeConProt(linkSiteList):
    k = 0
    while k < len(linkSiteList):
        if "CON_" in linkSiteList[k]:
            del linkSiteList[k]
        else:
            k += 1


def reportXLfile(openedXLfile, flTOwrt, linktype):
    f = openedXLfile; b = flTOwrt
    i = 2
    while i < len(f):
        print(i)
        lineList = f[i].split(",")
        linkPair = [lineList[1]]
        numSpec = lineList[-1].strip()
        m = i + 1
        if "SameSet" in f[m] or "SubSet" in f[m]:
            while m < len(f):
                if f[m].split(",")[0] == "":
                    break
                else:
                    linkPair.append(f[m].split(",")[1])
                    m += 1
        removeConProt(linkPair)
        linksite = "/".join(linkPair)
        bestSVM, bestEvale, pep, spec1, spec2, endPos = continueLines(m, f, numSpec) 
        if linksite != "":               
            wList = [linksite, numSpec, bestEvale, bestSVM, pep, linktype, spec1, spec2]
            b.write(",".join(wList) + "\n")
        i = endPos


def main():
    for fl in os.listdir("./"):
        if fl.endswith("filtered_cross-linked_sites.csv"):
            xlfl = open(fl, 'r').readlines()
        if fl.endswith("filtered_loop-linked_sites.csv"):
            lpfl = open(fl, 'r').readlines()
    b = open("report_chem.csv", 'w')
    b.write(",".join(["linkPair", "numSpec", "bestEvale", "bestSVM", "pep", "LinkType", "spec1", "spec2"])+"\n")
    reportXLfile(xlfl, b, "Cross-link")
    reportXLfile(lpfl, b, "Loop-link")
    b.close()


if __name__ == "__main__":
    main()