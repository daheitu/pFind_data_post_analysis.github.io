# coding = utf-8

def get_target():

    f = open(r"C:\Users\Yong Cao\Downloads\Belsom_Rappsilber_sulfoSDA.csv", 'r').readlines()

    b =open("report.csv", 'w')

    spec_dic = {}
    pep_dic = {}
    for line in f[1:]:
        linelist = line.split(",")
        raw = linelist[1]
        xlink = linelist[8:12]
        print(linelist[8:12])
        if raw not in spec_dic:
            spec_dic[raw] = 1
        else:
            spec_dic[raw] += 1
        
        if raw not in pep_dic:
            pep_dic[raw] = [xlink]
        else:
            if xlink not in pep_dic[raw]:
                pep_dic[raw].append(xlink)
        
    print(spec_dic)

    raw_list = sorted(list(spec_dic.keys()), key = lambda x: int(x.split('_')[-1]))

    for raw in raw_list:
        b.write(",".join([raw, str(len(pep_dic[raw])), str(spec_dic[raw])])+"\n")
        # print(raw, len(pep_dic[raw]))

    a = sorted(pep_dic.items(), key = lambda x:len(x[1]), reverse = True)
    target = [x[0] for x in a[:10]]
    return target