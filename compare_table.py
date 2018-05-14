import os

os.chdir(r"C:\Users\Yong\Desktop\results\combine")

f = open("report.txt", 'r').readlines()
b = open("compare.txt", 'w')
for line in f:
    line_list = line.strip().split("\t")
    count_list = [line_list[1], line_list[5]]
    for i in range(len(count_list)):
        if count_list[i] == "":
            count_list[i] = "0"
        else:
            continue
    
    new_list = []
    if int(count_list[0]) > int(count_list[1]):
        new_list = line_list[:5]
        new_list.insert(2, "")
        new_list.append(line_list[-1])
    else:
        new_list = line_list[5:]
        new_list.insert(0, line_list[0])
    print(new_list)
    b.write("\t".join(new_list))
    b.write("\n")
b.close()