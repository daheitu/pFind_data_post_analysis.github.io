fasta1 = r"K:\LiuFan_FAIMS_AC\DSSO\HEK_DSSO_BioRep1.fasta"
fasta2 = r"K:\LiuFan_FAIMS_AC\DSSO\HEK_DSSO_BioRep2.fasta"



def readfasta_info(fatsa):
    fas_dic = {}
    f = open(fatsa).readlines()
    i = 0 
    while i < len(f):
        if f[i][0] != ">":
            i += 1
        else:
            pro = f[i].rstrip("\n")
            seq = ""
            p = i + 1
            while p < len(f):
                if f[p][0] == ">":
                    break
                else:
                    seq += f[p]
                    p += 1
            fas_dic[pro] = seq
            i = p
    return fas_dic


fas1_dic = readfasta_info(fasta1)
fas2_dic = readfasta_info(fasta2)

for name in fas2_dic:
    if name not in fas1_dic:
        print(name)
        fas1_dic[name] = fas2_dic[name]

b = open("HEK_293_DSSO.fasta", 'w')
for name in fas1_dic:
    b.write(name+"\n")
    b.write(fas1_dic[name])
b.close()