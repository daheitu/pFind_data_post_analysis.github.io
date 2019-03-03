import os

os.chdir(r"G:\DSSO_1008\CID_BASED_Methods\ETD")

def generate_ion_mass_range(num, tol_ppm):
    deta = num * tol_ppm / 1000000
    return num - deta, num + deta


def compareTwoMZ(mz1, mz2, tol_ppm):
    lowMZ1,upMZ1 = generate_ion_mass_range(mz1, tol_ppm)
    if mz2 > lowMZ1 and mz2 < upMZ1:
        return True
    else:
        return False


def readETDinfo(path):
    etd = open(path, 'r').readlines()
    info_dic = {}
    i = 0
    while i < len(etd):
        if etd[i].strip() == "BEGIN IONS":
            title = etd[i+1].strip().split("=")[1]
            scanNum = int(title.split('.')[1])
            charge = int("".join([i for i in etd[i+2] if i.isdigit()]))
            mz = etd[i+4].strip().split("=")[1]
            info_dic[scanNum] = [charge, mz]
        else:
            print('wrong')
        p = i+5
        while p < len(etd):
            if etd[p].strip() == "END IONS":
                break
            else:
                pass
            p += 1
        i = p + 1
    return info_dic

# print(info_dic)

def readCIDinfo(path):
    cid = open(path, 'r').readlines()
    info_dic = {}
    i = 0
    while i < len(cid):
        if cid[i].strip() == "BEGIN IONS":
            title = cid[i+1].strip().split("=")[1]
            scanNum = int(title.split('.')[1])
            charge = int("".join([i for i in cid[i+2] if i.isdigit()]))
            mz = cid[i+4].strip().split("=")[1]
            info_dic[scanNum] = [charge, mz]
        else:
            print('wrong')
        p = i+5
        while p < len(cid):
            if cid[p].strip() == "END IONS":
                break
            else:
                pass
            p += 1
        i = p + 1
    return info_dic


def main():
    repList=[]
    delta_num = []
    etdDic = readETDinfo("DSSO_CID_ETD_MS2_RE23_MIPS_T1_ETDFT.mgf")
    cidDic = readETDinfo("DSSO_CID_ETD_MS2_RE23_MIPS_T1_CIDFT.mgf")
    for num in etdDic:
        etdCharge, etdMZ = etdDic[num]
        for k in range(1, 15):
            if num + k not in cidDic:
                pass
            else:
                cidCharge, cidMZ = cidDic[num+k]
                if etdCharge == cidCharge and compareTwoMZ(float(etdMZ), float(cidMZ), 10):
                    pair = [num, num+k]
                    repList.append(pair)
                    delta_num.append(k)
                    break
                else:
                    pass
    delta_num.sort()
    print(repList)
    print(delta_num)


if __name__ == "__main__":
    main()
