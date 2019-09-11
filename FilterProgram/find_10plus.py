f = open('./min_combin_4.txt', 'r').readlines()
b = open("10plus.txt", 'w')
for line in f:
    if float(line.strip().split("\t")[-1]) >10:
        b.write(line)
b.close()