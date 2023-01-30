# coding = utf-8

import os

# wkdir = r"K:\20200801"


def get_link_sites(rpfl_path):
    sites_list = []
    f = open(rpfl_path).readlines()
    for line in f[1:]:
        lineList = line.split(",")
        if lineList[0].isdigit():
            break
        else:
            sites_list.append(lineList[0])
    return sites_list


def stac_final(wkdir, rep_dic):
    # rep_dic = {}
    for root, dirs, fls in os.walk(wkdir):
        if root.endswith("reports"):
            linker, out_put = root.split("\\")[-3:-1]
            if linker == "DSSO":
                if out_put == "output_plink_dsso":
                    for fl in fls:
                        if fl.endswith("v5.csv"):
                            rep_path = os.path.join(root, fl)
                            if linker not in rep_dic:
                                rep_dic[linker] = get_link_sites(rep_path)
                            else:
                                rep_dic[linker].extend(get_link_sites(rep_path))
            else: 
                if out_put == "output":
                    # if # linker == "DSS":
                    for fl in fls:
                        if fl.endswith("v5.csv"):
                            rep_path = os.path.join(root, fl)
                            if linker not in rep_dic:
                                    rep_dic[linker] = get_link_sites(rep_path)
                            else:
                                site_list = get_link_sites(rep_path)
                                for site in site_list:
                                    if site not in rep_dic[linker]:
                                        rep_dic[linker].append(site)
                                # rep_dic[linker].extend(get_link_sites(rep_path))



def find_site_v5_file(flpath, site):
    f = open(flpath).readlines()
    for line in f:
        lineList = line.split(",")
        if lineList[0].isdigit():
            return False, None
        elif lineList[0] == site:
            info_list = []
            for i in range(6, len(lineList), 3):
                info_list.append(lineList[i])
            return True, info_list
        else:
            continue
    return False, None


def isright_folder(root, lk):
    if root.endswith("reports"):
        # print(root, lk)
        linker, out_put = root.split("\\")[-3:-1]
        if linker == lk:
            if linker == "DSSO":
                if out_put == "output_plink_dsso":
                    return True
                else:
                    return False
        
            else:
                if out_put == "output":
                    return True
                else:
                    return False
                
        else:
            return False
    else:
        return False


def find_site(wkpath, lk, site):
    for root, dirs, fls in os.walk(wkpath):
        # print(root)
        if isright_folder(root, lk):
            # print("%s is right" % root)
            for fl in fls:
                if fl.endswith("v5.csv"):
                    rep_path = os.path.join(root, fl)
                    isfind, info = find_site_v5_file(rep_path, site)
                    if isfind:
                        return True, info
                    else:
                        continue
    return False, None


# print(rep_dic)
def main_flow():
    rep_dic = {}
    path1 = r"G:\msData\20200419"
    path2 = r"G:\msData\20200529"
    stac_final(path1, rep_dic)
    stac_final(path2, rep_dic)
    over_set = set(rep_dic["DSS"])
    for lk in rep_dic:
        over_set = over_set | set(rep_dic[lk])
    print(len(over_set))
    b = open("merge_4_all.csv", 'w')
    for site in list(over_set)[:]:
        wlist = [site]
        for lk in ["DSS", "BSMEG", "DSS_C", "DSSO"]:
            for path in [path1, path2]:
                isfind, info = find_site(path, lk, site)
                if isfind:
                    wlist.extend(info)
                    break
                else:
                    continue
            if not isfind:
                wlist.extend([""]*4)
        # print(wlist)
        b.write(",".join(wlist)+"\n")
    b.close()
    



    # b = open("total_summary_4.csv", 'w')
    # # print(rep_dic)
    # row = len(rep_dic["DSS"])
    # for i in range(row):
    #     # if i > len(rep_dic[])
    #     wlist = []
    #     for linker in ["DSS", "BSMEG", "DSS_C", "DSSO"]:#"
    #         if i >= len(rep_dic[linker]):
    #             wlist.append("")
    #         else:
    #             wlist.append(rep_dic[linker][i])
    #     b.write(",".join(wlist)+"\n")
    # b.close()


main_flow()