b = open("stepHCD", 'w')
for i in range(21, 25):
    for j in range(30,34):
        a = (i + j) 
        if a % 2 == 0:
            b.write("[%d %d %d]" %(i, int(a/2), j) + ",")
b.close()