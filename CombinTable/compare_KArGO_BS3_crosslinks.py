import os
import re

os.chdir(r"C:\Users\Yong\Desktop\UTPA")
cut_off_min = 8
cut_off_max = 20

def type_judgement(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find("-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:n].strip()
    if protein1 != protein2:
        link_type = "Inter"
    else:
        if abs(int(position1) - int(position2)) < 5:
            link_type = "Inter"
        else:
            link_type = "Intra"
    return [protein1, position1, protein2, position2, link_type]


def compare_KK_KR_links(linkpair1, linkpair2):
    linktype1 = type_judgement(linkpair1)[-1]
    linktype2 = type_judgement(linkpair2)[-1]
    pos_list_1 = [int(type_judgement(linkpair1)[1]), int(type_judgement(linkpair1)[3])]
    pos_list_2 = [int(type_judgement(linkpair2)[1]), int(type_judgement(linkpair2)[3])]
    pro_list_1 = [type_judgement(linkpair1)[0], type_judgement(linkpair1)[2]]
    pro_list_2 = [type_judgement(linkpair2)[0], type_judgement(linkpair2)[2]]
    if linktype1 != linktype2:
        return False
    else:
        if linktype1 == "Intra" and pro_list_1 == pro_list_2:
            delta_list = [abs(pos_list_2[0]-pos_list_1[0]), abs(pos_list_2[0]-pos_list_1[1]), abs(pos_list_2[1]-pos_list_1[0]), abs(pos_list_2[1]-pos_list_1[1])]
            # print(delta_list)
            if min([delta_list[0], delta_list[3]]) <= min(delta_list[1], delta_list[2]):
                if min([delta_list[0], delta_list[3]]) < cut_off_min and max([delta_list[0], delta_list[3]]) < cut_off_max:
                    return True
                else:
                    return False
            else:
                if min([delta_list[1], delta_list[2]]) < cut_off_min and max([delta_list[1], delta_list[2]]) < cut_off_max:
                    return True
                else:
                    return False
        else:
            if pro_list_2[0] == pro_list_1[0] and pro_list_2[1] == pro_list_1[1]:
                if abs(pos_list_1[0] - pos_list_2[0]) < cut_off_max and abs(pos_list_1[1] - pos_list_2[1]) < cut_off_max:
                    return True
                else:
                    return False
            elif pro_list_2[0] == pro_list_1[1] and pro_list_2[1] == pro_list_1[0]:
                if abs(pos_list_1[1] - pos_list_2[0]) < cut_off_max and abs(pos_list_1[0] - pos_list_2[1]) < cut_off_max:
                    return True
                else:
                    return False
            else:
                return False
print(compare_KK_KR_links("Utp4 (313)-Utp5 (41)", "Utp4(313)-Utp5 (58)"))

def main():
    f1 = open("BS3_Inter_intra.txt").readlines()
    f2 = open("KArGO_inter_intra.txt").readlines()
    rep = open("report.txt", 'w')
    bs3_dic = {}
    for line in f1:
        line_list = line.strip().split("\t")
        linkpair = line_list[0]
        if linkpair not in bs3_dic:
            bs3_dic[linkpair] = line_list
        else:
            continue
    
    KArGO_dic = {}
    for line in f2:
        line_list = line.strip().split("\t")
        links = line_list[0]
        if links not in KArGO_dic:
            KArGO_dic[links] = line_list
        else:
            continue
    
    for linkpair in bs3_dic:
        for links in list(KArGO_dic.keys()):
            if compare_KK_KR_links(linkpair, links) is True:
                # print(linkpair, links, compare_KK_KR_links(linkpair, links))
                #print(bs3_dic[linkpair], KArGO_dic[links])
                bs3_dic[linkpair] += KArGO_dic[links]
                #print(bs3_dic[linkpair])
            else:
                continue

    for linkpair in bs3_dic:
        rep.write("\t".join(bs3_dic[linkpair]))
        rep.write("\n")
    
    rep.close()

if __name__ == "__main__":
    main()
    