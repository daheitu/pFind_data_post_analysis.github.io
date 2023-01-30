# coding = utf-8
import os
os.chdir(r"E:\Script\test_SGC_20190512") # pFind.protein 所在文件夹
modify = "Smo_1"  # important! enter the modification type you want search


def filteModpFind(openedFl, idxPepStar, idxPepEnd):
    repDic = {}
    finalDic = {}
    tab = openedFl
    for line in tab[idxPepStar:idxPepEnd +1]:
        line_list = line.rstrip("\n").split("\t")
        if modify not in line_list[8]:
            continue
        else:
            site_info = []
            pepSeq = line_list[4]
            modi_list = line_list[8].split(";")[:-1]
            protList = line_list[10].split("/" )[:-1]
            pepStarPosList = []
            for pepPos in line_list[11].split("/" )[:-1]:
                starPos = int(pepPos.split(",")[0])
                pepStarPosList.append(starPos)
            for mod in modi_list:
                mod_type = mod.split(",")[1]
                modSite = int(mod.split(",")[0])
                if modify == mod_type:  # modi name
                    modSite = int(mod.split(",")[0])
                    for i in range(len(protList)):
                        prot = protList[i]
                        pepStarPos = pepStarPosList[i]
                        modSiteInPro = str(pepStarPos + modSite)
                        site_info.append(prot + "[" + modSiteInPro + "]")
        
            if site_info:
                SITE_IFOR = ";".join(site_info)
                peptide_ID = line_list[2]
                peptide_seq = line_list[3]
                final_score = float(line_list[7])
                peptide_mod = line_list[8]
                peptide_spec_num = int(line_list[-1])
                protein_name = line_list[10]
                spectra_title = line_list[-3]
                repDic[peptide_ID] = ["", peptide_ID, peptide_seq, final_score,\
                    peptide_mod, peptide_spec_num, protein_name, \
                    spectra_title, SITE_IFOR]
    if repDic == {}:
        return None
    else:
        for order in repDic:
            SITE_IFOR = repDic[order][-1]
            site_info = SITE_IFOR.split(";")
            peptide_spec_num = repDic[order][5]
            final_score = repDic[order][3]
            for site in site_info:
                if site not in finalDic:
                    finalDic[site] = ["", "", site, peptide_spec_num, final_score]
                else:
                    finalDic[site][-2] += peptide_spec_num
                    finalDic[site][-1] = min(finalDic[site][-1], final_score)
        
        return repDic, finalDic


def wirteDicToFl(wtDic, flToW):
    for key in wtDic:
        wList = [str(ele) for ele in wtDic[key]]
        flToW.write("\t".join(wList) + "\n")


def main():
    f = open("pFind.protein", 'r').readlines()
    b = open("report.txt", 'w')
    b.write("\t".join([
        "ID", "Sequence", "Protein", "Modification", "Modi_site",
        "Spectra Number", "Final score", "Spectra title"
    ]))
    b.write("\n")
    i = 2
    while i < len(f):
        if f[i].rstrip("\n").split("\t")[0].isdigit():
            if "REV_" in f[i]:
                m = i + 1
                while m < len(f):
                    if f[m].rstrip("\n").split("\t")[0].isdigit():
                        break
                    else:
                        m += 1
                i = m
            else:
                idxStar = i
                p = i + 1
                while p < len(f):
                    #lineList = f[p].rstrip("\n").split("\t")
                    if f[p].rstrip("\n").split("\t")[0] == "" and \
                         f[p].rstrip("\n").split("\t")[2] == "1":
                        idxPepStar = p
                        p += 1
                    elif f[p].rstrip("\n").split("\t")[0].isdigit():
                        idxPepEnd = p - 1
                        break
                    elif f[p][:4] == "----":
                        idxPepEnd = p - 1
                        break
                    else:
                        p += 1

                print(idxStar, idxPepStar, idxPepEnd)
                if filteModpFind(f, idxPepStar, idxPepEnd) != None:
                    infoDic, statsDic = filteModpFind(f, idxPepStar, idxPepEnd)
                    for line in f[idxStar: idxPepStar]:
                        if "REV_" not in line:
                            b.write(line)
                    wirteDicToFl(infoDic, b)
                    wirteDicToFl(statsDic, b)
                    
                i = p
    
        elif f[i][:4] == "----":
            break
    b.close()


if __name__ == "__main__":
    main()
