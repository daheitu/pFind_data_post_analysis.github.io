import re


def getSiteInfo(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]
    return protein1, position1, protein2, position2


def judegGUI(protein, p1F):
    if protein[0] == p1F[0] and protein.endswith(p1F[1:]):
        return True
    else:
        return False


def judgesiteGUI(protein1, p1F, p2F):
    if judegGUI(protein1, p1F):
        return "p1F"
    elif judegGUI(protein1, p2F):
        return "p2F"
    else:
        return False


def calsyntheticNum(fileName, p1F, p2F):
    f = open(fileName, 'r').readlines()
    interCrossmid = 0
    interCrossmid_spec = []
    interCrossNter = 0
    interCrossNter_spec = []
    intraCrossmid = 0
    intraCrossmid_spec = []
    intraCrossNter = 0
    intraCrossNter_spec = []
    for line in f[1:]:
        lineList = line.split(",")
        linkPair = lineList[0]
        spec = int(lineList[1])
        pro1, site1, pro2, site2 = getSiteInfo(linkPair)
        protype1 = judgesiteGUI(pro1, p1F, p2F)
        protype2 = judgesiteGUI(pro2, p1F, p2F)
        if protype1 and protype2:
            if "p1F" in [protype1, protype2] and "p2F" in [protype1, protype2]:
                if site1 != "1" and site2 != "1":
                    interCrossmid += 1
                    interCrossmid_spec.append(spec)
                else:
                    interCrossNter += 1
                    interCrossNter_spec.append(spec)
            else:
                if site1 != "1" and site2 != "1":
                    intraCrossmid += 1
                    intraCrossmid_spec.append(spec)
                else:
                    intraCrossNter += 1
                    intraCrossNter_spec.append(spec)
                       
    print("interCrossmid is %d\ninterCrossNter is %d\nintraCrossmid is %d\nintraCrossNter is %d" % (interCrossmid, interCrossNter, intraCrossmid, interCrossNter))
    print("The spec of interCrossmid is %d" % sum(interCrossmid_spec))
    print("The spec of interCrossNter is %d" % sum(interCrossNter_spec))
    print("The spec of intraCrossmid is %d" % sum(intraCrossmid_spec))
    print("The spec of intraCrossNter is %d" % sum(intraCrossNter_spec))


def calsyntheticNumforFix(fileName, p1F):
    f = open(fileName, 'r').readlines()
    interCrossmid = 0
    interCrossmid_spec = []
    interCrossNter = 0
    interCrossNter_spec = []
    intraCrossmid = 0
    intraCrossmid_spec = []
    intraCrossNter = 0
    intraCrossNter_spec = []
    for line in f[1:]:
        lineList = line.split(",")
        linkPair = lineList[0]
        spec = int(lineList[1])
        pro1, site1, pro2, site2 = getSiteInfo(linkPair)
        protype1 = judegGUI(pro1, p1F)
        protype2 = judegGUI(pro2, p1F)
        if protype1 and protype2:
            if site1 != "1" and site2 != "1":
                intraCrossmid += 1
                intraCrossmid_spec.append(spec)
            else:
                intraCrossNter += 1
                intraCrossNter_spec.append(spec)
        elif protype1 == False and protype2 == False:
            pass
        else:
            if "Fix" in pro1 or "Fix" in pro2:
                if site1 != "1" and site2 != "1":
                    interCrossmid += 1
                    interCrossmid_spec.append(spec)
                else:
                    interCrossNter += 1
                    interCrossNter_spec.append(spec)
                       
    print("interCrossmid is %d\ninterCrossNter is %d\nintraCrossmid is %d\nintraCrossNter is %d" % (interCrossmid, interCrossNter, intraCrossmid, interCrossNter))
    print("The spec of interCrossmid is %d" % sum(interCrossmid_spec))
    print("The spec of interCrossNter is %d" % sum(interCrossNter_spec))
    print("The spec of intraCrossmid is %d" % sum(intraCrossmid_spec))