from Format_FastaFile import formatFast
f = open(r"C:\Users\Yong Cao\Documents\WeChat Files\yunzaitianya2010\FileStorage\File\2019-09\NPC_9seq.fasta").readlines()
b = open("NPL.fasta", 'w')
i = 0
while i < len(f):
    if f[i][0] != ">":
        i += 1
    else:
        b.write(f[i])
        p = i+ 1
        while p < len(f):
            if f[p][0] == ">":
                break
            else:
                a = f[p].rstrip().upper()
                b.write(a + "\n")
                p +=1
        i = p
b.close()
formatFast("./X.laevis_Xenbase.fasta")