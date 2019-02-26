import re

database_path = r"C:\Users\Yong\Desktop\atg16_database\human_unitprot_contaminant.fasta"

mf2 = "H[GAVLISTKDNEQPWFYM]{5,9}R"


fasta = open(database_path, 'r').readlines()
b = open("report.txt", 'w')
i = 0
while i < len(fasta) and fasta[i][0] == '>':
    # print(i)
    if fasta[i][1:5] == "CON_":
        for m in range(len(fasta[i])):           
            if fasta[i][m] == "|":  # in [" ", "|", "\n", "\t"]:
                line_Pro_name = fasta[i][1:m]
                break
            else:
                continue
    else:
        for m in range(len(fasta[i])):           
            if fasta[i][m] in [" ", "\n", "\t"]:
                line_Pro_name = fasta[i][1:m]
                break
            else:
                continue
    p = i + 1

    pro_seq = ""
    while p < len(fasta):
        if fasta[p][0] == '>':
            break
        else:
            pro_seq += fasta[p].strip()
        p += 1
    pos = re.findall(mf2, pro_seq)
    i = p
    if pos:
        for pep in pos:
            if pep.count("S") == 1:
                print(pep, line_Pro_name)
            else:
                continue
    else:
        continue
      
b.close()