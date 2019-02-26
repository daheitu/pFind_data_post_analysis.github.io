import os
import re
os.chdir(
    r"C:\Users\Yong\Documents\pLink\pLink_task_2018.05.22.15.46.42_BSA_BUFFER_SCREEN\reports"
)
file_list = os.listdir(os.getcwd())

protein_list = ["Nop10 ", "Cbf5 ", "Gar1 ", "Nhp2 ", "lysozyme", "BSA", "GST", "lactoferrin", "PUD1", "PUD2", "Aldolase"]

raw_list = ["Tris(BORATE)", "Tris_50mm", "borae", "hepes(BORATE)", "hepes_50"]

def judge_sites(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]
    if protein1 in protein_list and protein2 in protein_list:
        return True
    else:
        return False



for file in file_list:
    if "cross-linked_peptides" in file:
        print(file)
        f = open(file, 'r').readlines()
        pep_dic = {}
        pep_list_pos = []
        pep_spec_list = []
        n = 2
        while n < len(f):
            line = f[n]
            line_list = line.rstrip("\n").split(",")
            if line_list[0].isdigit():
                pep = line_list[1]
                link_site = line_list[4]
                p = n+1
            else:
                print(n)
            sub_pep = {}
            while p < len(f) and f[p].rstrip("\n").split(",")[0] == "":
                for word in raw_list:
                    if word in f[p]:
                        if word not in sub_pep:
                            sub_pep[word] = 1
                        else:
                            sub_pep[word] += 1
                        break
                    else:
                        continue
                p += 1
            
            print(sub_pep)
            for word in sub_pep:
                if sub_pep[word] > 1:
                    if word not in pep_dic:
                        pep_dic[word] = 1
                    else:
                        pep_dic[word] += 1
                else:
                    continue
        
            n = p
        print(pep_dic)
        """
        b = open("pep_report.txt", "w")
        num_0dot1 = 0
        num_0dot2 = 0
        for key in pep_dic:         
            if judge_sites(pep_dic[key][1][:-1]):
                if int(pep_dic[key][2]) > 1 and int(pep_dic[key][3]) > 1:
                    b.write("\t".join(pep_dic[key]))
                    b.write("\n")
                    num_0dot1 += 1
                    num_0dot2 += 1
                elif int(pep_dic[key][2]) > 1 and int(pep_dic[key][3]) < 2:
                    b.write("\t".join(pep_dic[key]))
                    b.write("\n")
                    num_0dot1 += 1
                elif int(pep_dic[key][2]) < 2 and int(pep_dic[key][3]) > 1:
                    b.write("\t".join(pep_dic[key]))
                    b.write("\n")
                    num_0dot2 += 1
            else:
                continue
        print(num_0dot1, num_0dot2)
"""