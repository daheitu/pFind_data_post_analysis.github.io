import os

os.chdir(r"D:\workspace\pFindTask60_ALPK1_NC\result")
f = open('pFind.protein', 'r').readlines()
count_N = 0
count_C = 0
count_midle = 0
for line in f[175:2861]: #2861, 1639
    line_list = line.strip().split("\t") 
    if int(line_list[9].split(',')[0]) < 500:
        count_N += int(line_list[-1])
    elif int(line_list[9].split(',')[0]) > 900:
        count_C += int(line_list[-1])
    else:
        count_midle += int(line_list[-1])

print(count_N, count_C, count_midle)

