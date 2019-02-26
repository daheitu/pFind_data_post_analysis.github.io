f = open(r'./report.csv', 'r').readlines()
name_list = []
for line in f[1:]:
    name_list.append(line.split(",")[0])

print(name_list)