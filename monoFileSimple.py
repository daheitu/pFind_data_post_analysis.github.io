import os

os.chdir(r"C:\Users\Yong\Documents\pLink\pLink_task_2019.03.13.14.23.50\reports")

def getMonoSiteFL(folderpath):
    flList = os.listdir(folderpath)
    for fl in flList:
        if fl[-21:] == "mono-linked_sites.csv":
            return fl
        else:
            continue
    if fl == flList[-1]:
        return None

def main():
    monofl = getMonoSiteFL(os.getcwd())
    f = open(monofl, 'r').readlines()
    b = open("monoLinkSiteReport.csv", 'w')
    i = 0
    while i < len(f):
        if not f[i][0].isdigit():
            i += 1
        else:
            lineList = f[i].strip().split(",")
            order = lineList[0]
            site = lineList[1]
            numOfPep = lineList[2]
            numOfSep = lineList[3]
            p = i + 1
            while p < len(f):
                if f[p].split(',')[0] != "":
                    p += 1
                else:
                    pep = f[p].split(",")[5]
                    bestSVM = f[p].split(",")[9]
                    break
            wlist = [order, site, numOfSep, numOfPep, pep, bestSVM]
            b.write(",".join(wlist) + "\n")
            i += 2


if __name__ == "__main__":
    main()