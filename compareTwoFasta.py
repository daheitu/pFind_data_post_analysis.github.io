import os
os.chdir(r"F:\MS_DATA_STORAGE\opti_DSSO_JPR2019")
pathK12 = r"D:\fasta database\uniprot_Ecoli_K12.fasta"
pathRibo = r"F:\MS_DATA_STORAGE\opti_DSSO_JPR2019\Ribosome_for_XL_all_contained_Proteins.fasta"



def getProteinNameSet(filepath):
    with open(filepath) as file:
        f = file.readlines()
    K12List = []
    for line in f:
        if line[0] == ">":
            name = line[1:line.find("|", 4)]
            if name not in K12List:
                K12List.append(name)
    return set(K12List)

k12set = getProteinNameSet(pathK12)
riboset = getProteinNameSet(pathRibo)

if riboset & k12set == riboset:
    print("OK")
else:
    print(riboset & k12set ^ riboset)

b = open("ribo_madeup.fasta", 'w')
fribo = open(pathRibo, 'r').readlines()
for line in fribo:
    b.write(line)

fk12 = open(pathK12, 'r').readlines()

i = 0
while i < len(fk12):
    if fk12[i][0] != ">":
        i += 1
    else:
        name = fk12[i][1:fk12[i].find("|", 4)]
        #print(name)
        if name in list(riboset):
            i += 1
        else:
            b.write(fk12[i])
            p = i+1
            while p < len(fk12):
                if fk12[p][0] == ">":
                    break
                else:
                    b.write(fk12[p])
                    p += 1
            i = p

b.close()