"""
f = open("./nmeth.3283-S7.csv", 'r').readlines()
print(len(f[4].rstrip("\n").split(",")))
b = open("report.csv", 'w')
name_list = []
for line in f:
    line_list = line.rstrip("\n").split(",")
    # print(len(line_list))
    name = line_list[0]
    name_list.append(name)
# print(name_list)
name_dic = {}
for name in name_list:
    if name not in name_dic:
        name_dic[name] = name_list.count(name)
    else:
        pass

for name in name_dic:
    if name_dic[name] > 5:
        b.write(name+ "\n")
b.close()
"""
b  = open("CFH-DISULFIDE.csv", 'w')
site_list = []
pCFH = [['21', '66'], ['52', '80'], ['85', '129'], ['114', '141'], ['146', '192'], ['178', '205'], ['210', '251'], ['237', '262'], ['267', '309'], ['294', '320'], ['325', '374'], ['357', '385'], ['389', '431'], ['416', '442'], ['448', '494'], ['477', '505'], ['509', '553'], ['536', '564'], ['569', '611'], ['597', '623'], ['630', '673'], ['659', '684'], ['691', '733'], ['719', '744'], ['753', '792'], ['781', '803'], ['811', '853'], ['839', '864'], ['870', '915'], ['901', '926'], ['931', '973'], ['959', '984'], ['989', '1032'], ['1018', '1043'], ['1048', '1091'], ['1077', '1102'], ['1109', '1152'], ['1138', '1163'], ['1167', '1218'], ['1201', '1228']]
for pair in pCFH:
    for site in pair:
        site_list.append(int(site))
#print(site_list)
site_list.sort()
# print(site_list)
for i in range(len(site_list)-1):
    if site_list[i+1] - site_list[i] < 6:
        if [str(site_list[i]), str(site_list[i+1])] not in pCFH:
            print(site_list[i+1], site_list[i])
            write_list = [str(site_list[i]), str(site_list[i+1])]
            for pair in pCFH:
                if str(site_list[i]) in pair:
                    print(pair)
                    write_list.append(str(pair))
                if str(site_list[i+1]) in pair:
                    print(pair)
                    write_list.append(str(pair))
        b.write(",".join(write_list)+"\n")
b.close()
