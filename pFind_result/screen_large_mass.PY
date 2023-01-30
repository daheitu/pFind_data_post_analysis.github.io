f = open(r"E:\pFindStudio\pFind316\bin\modification.ini", 'r').readlines()
b = open("unimod_E.txt", 'w')
i = 1
n = 0
while i in range(len(f)):
    if "[S]" in f[i]:
        if float(f[i+1].split(" ")[2]) > 100:
            n += 1
            b.write(f[i+1])
            print(f[i+1])
    i += 2
        
print(n)    
b.close()