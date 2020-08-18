import os
import re

wk_dir = r"K:\20200801"
def judgeHomo(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find(")-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m]
    protein2 = linked_site[p + 2:n]
    if protein1 == protein2:
        return protein1.strip()
    else:
        return False


def get_protein_list(wk_dir):
    protein_list = []
    for root, dirs, fls in os.walk(wk_dir):
        for fl in fls:
            if fl.endswith("cross-linked_sites.csv"):
                print(os.path.join(root, fl))
                f = open(os.path.join(root, fl)).readlines()
                for line in f[2:]:
                    linelist = line.split(',')
                    if len(linelist) == 4:
                        site = linelist[1]
                        protein = judgeHomo(site)
                        if protein:
                            if protein not in protein_list:
                                protein_list.append(protein)
    return protein_list
# print(len(protein_list))
protein_list = get_protein_list(wk_dir)
def write_fasta(protein_list):
    print(len(protein_list))
    fasta_path = r"D:\FastaDatabase\Ecoli-uniprot-mg1655-20160918.fasta"
    f = open(fasta_path).readlines()
    b = open("small.fasta", 'w')
    i = 0 
    while i < len(f):
        # print(f[i])
        if f[i][0] != ">":
            i += 1
        else:
            name = f[i].split(" ")[0][1:].strip()
            if name not in protein_list:
                i += 1
            else:
                # print(i)
                b.write(f[i])
                i += 1
                while i < len(f):
                    if f[i].startswith(">"):
                        break
                    else:
                        b.write(f[i])
                        i += 1
                
    b.close()
write_fasta(protein_list)
f = open('./small.fasta').readlines()
small_list = []
for line in f:
    if line.startswith(">"):
        name = line.split(" ")[0][1:].strip()
        small_list.append(name)
for p in protein_list:
    if p not in small_list:
        print(p)