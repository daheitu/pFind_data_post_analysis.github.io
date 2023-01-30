
import os
os.chdir(r"E:\workspace\pFindTask114_YANZ_trans_id\result")
list_file = "./transcolon_BS3_proteinID.csv"
pFind_results = "./pFind.protein"

topN = 3000

def find_group_leader(group_name, f, line_index):
    want_lines = []
    i = 0
    while i < len(f):
        line_list = f[i].split("\t")
        if line_list[0].isdigit() and line_list[1] == group_name:
            line_mod_list = [str(line_index)]
            line_mod_list.extend(line_list[1:])
            want_lines.append("\t".join(line_mod_list))
            i += 1
            while i < len(f):
                if f[i].split("\t")[0].isdigit():
                    break
                else:
                    want_lines.append(f[i])
                    i += 1
        else:
            i += 1
    return want_lines


# pro_file = open(pFind_results, 'r').readlines()
# f = open(list_file, 'r').readlines()
# c = open("./pFind_protein_target_protein.txt", 'w')
# c.write(pro_file[0])
# c.write(pro_file[1])
# for line in f[1:topN+1]:
#     line_index = f.index(line)
#     line_list = line.split(",")
#     group_name = line_list[1]
#     group_leader = find_group_leader(group_name, pro_file, line_index)
#     print(group_leader)
#     c.writelines(group_leader)

# c.close()

def find_group_peptides_num(group_name):
    f = open(pFind_results).readlines()
    for line in f:
        if line[0].isdigit():
            line_list = line.split("\t")
            # print(line_list)
            protein = line_list[1]
            if protein == group_name:
                return line_list[5]
    return "0"


tgt = open(list_file).readlines()
rep = open("protein_id.csv", 'w')
rep.writelines(tgt[0])
for line in tgt[1:]:
    line_list = line.split(',')
    prot = line_list[1]
    pep_num = find_group_peptides_num(prot)
    line_list.insert(5, pep_num)
    rep.write(",".join(line_list))

rep.close()