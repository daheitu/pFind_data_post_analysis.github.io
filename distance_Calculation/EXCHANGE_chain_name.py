import os

os.chdir(r"E:\Rosetta Learning\prepack")
b = open("3u28_change.pdb", 'w')
f = open("3u28.pdb", 'r').readlines()

def writeInfo(chain, pdb, flToWT):
    for line in pdb:
        if line[:3] == "TER" and str(line[21]) == chain:
            flToWT.write(line)
            break
        elif line[:4] == "ATOM" and str(line[21]) == chain:
            flToWT.write(line)
        else:
            continue
def changeBToC(chainR, chainW, pdb, flToWT):
    for line in pdb:
        if line[:3] == "TER" and str(line[21]) == chainR:
            flToWT.write(line[:21]+chainW+line[22:])
            break
        elif line[:4] in ["ATOM", "ANIS"] and str(line[21]) == chainR:
            flToWT.write(line[:21]+chainW+line[22:])
        else:
            continue

writeInfo("A", f, b)
changeBToC("B", "C", f, b)
changeBToC("C", "B", f, b)

b.close()