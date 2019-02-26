import os
os.chdir(r"G:\DSSO_1008\structure")
chain = "B"
input_file = "3v03.pdb"

pdb = open(input_file, 'r').readlines()
b = open(input_file[:-4] + "_monomer.pdb", 'w')

for line in pdb:
    if str(line[21]) == chain:
        pass
    else:
        b.write(line)

b.close()