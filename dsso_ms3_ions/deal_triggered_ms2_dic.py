# coding = utf-8

from readMS3info import *

ms3InfoDic = readms3Info("./DSSO_CID_MS2_CID_MS3_T1.ms3")
numTotalMS3 = len(ms3InfoDic)
ms2TriggerMS3dic = getLinageInfo(ms3InfoDic)
totalms2, ms2TriggerMS3dic = readms2Info("./DSSO_CID_MS2_CID_MS3_T1.ms2", ms2TriggerMS3dic)
# writeMs2ToMs3Info(ms2TriggerMS3dic, "ms2_ms3_linage.txt")
mH = 1.00782
mS = 31.9720
mH2O = 18.0105
len_dic = {}

# for scan in ms2TriggerMS3dic:
#     # if ms2TriggerMS3dic[scan]:
#     # print(len(ms2TriggerMS3dic[scan]))
#     len_scan_list = len(ms2TriggerMS3dic[scan])
#     if len_scan_list not in len_dic:
#         len_dic[len_scan_list] = 1
#     else:
#         len_dic[len_scan_list] += 1

# print(len_dic)

# 根据一个ms3母离子的质量和ms2母离子质量判断是否为mono
# def isMono(ms2_info, ms3_ion_info):

# 判断两个离子是否成对
def isPaired(a_list, b_list, tolppm=30):
    mz_a, z_a = a_list[:2]
    mz_b, z_b = b_list[:2]
    if z_a != z_b:
        return False, 50
    else:
        if z_b != 0:
            min_mz = min(mz_a, mz_b)
            max_mz = max(mz_a, mz_b)
            ismatch, delta = compareTwoNum(min_mz+31.972/z_b, max_mz, tolppm)
            return ismatch, delta
        else:
            return False, 50


# 去冗余，将字典里成对离子中的大质量的删掉，保留小的
def remove_large_in_pair(all_ms3):
    non_re_info = []
    i = 0
    while i < len(all_ms3)-1:
        for j in range(i+1, len(all_ms3)):
            ismatch, delta = isPaired(all_ms3[i], all_ms3[j])
            # print(ismatch, delta)
            if ismatch:
                if all_ms3[i] < all_ms3[j]:
                    a = list(all_ms3[i])
                else:
                    a = list(all_ms3[j])
                a.append(delta)
                non_re_info.append(tuple(a))
                del all_ms3[i]
                del all_ms3[j-1]
                break
        if not ismatch:
            non_re_info.append(all_ms3.pop(i))
            # del all_ms3[i]

    if all_ms3 != []:
        non_re_info.append(all_ms3[-1])
    return non_re_info


def isFeatureIon(paired_ms3, ms2):
    mz1, z1 = paired_ms3[0][:2]
    mz2, z2 = paired_ms3[1][:2]
    mpre, zpre = ms2[1:3]
    mzpre_from_ms3 = ((mz1- mH)*z1 + (mz2-mH)*z2 + mS + mH2O + mH * zpre)/zpre
    ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
    # print(ismatch, delta)
    if ismatch:
        return True
    else:
        return False


# print(isFeatureIon([(522.7877, 2, 6023), (436.253, 2, 6025)], (6018, 655.68683, 3, 1965.04593)))
# print(remove_large_in_pair([(1524.7208, 3, 23394), (1535.3754, 3, 23395), (1535.0325, 3, 23396)]))

# '''
both_ms2_num = 0
for scan in ms2TriggerMS3dic:
    # if scan == 6018:
    ms2_info = ms2TriggerMS3dic[scan]
    if len(ms2_info) == 5:
        ms2 = ms2_info[0]
        all_ms3 = ms2_info[1:]
        # print(all_ms3)
        non_re_info = remove_large_in_pair(all_ms3)
        # print(non_re_info, ms2)
        paired_ms3 = [x for x in non_re_info if len(x) == 4]
        # print(paired_ms3)
        if len(paired_ms3) == 2:
            # print(paired_ms3)
            if isFeatureIon(paired_ms3, ms2):
                print(scan)
                print(paired_ms3, ms2)
                both_ms2_num += 1

print(both_ms2_num/len(ms2TriggerMS3dic))





# '''
        



    
            




