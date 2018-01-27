import os
os.chdir(r"C:\Users\Yong\Documents\pLink\search_task_2018.01.18.08.36.24_BSA_ARGO2_4RAW\reports")
file_list = os.listdir(os.getcwd())
for file in file_list:
    if "cross-linked_peptides" in file:
        print(file)
        f = open(file, 'r').readlines()
        n=0
        for line in f:
            if line.rstrip("\n").split(",")[0].isdigit():
                if float(f[f.index(line)+1].strip().split(",")[6]) < 0.5:
                    n +=1
        print(n)
