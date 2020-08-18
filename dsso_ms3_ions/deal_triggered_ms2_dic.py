# coding = utf-8

from readMS3info import readms3Info, getLinageInfo, readms2Info, compareTwoNum

<<<<<<< HEAD
wk_dir = r"F:\data_from_paper\maxlinker_MCP\JinLiang10366422"
os.chdir(wk_dir)

# for fl in os.listdir(wk_dir):
#     if fl.endswith(".ms2"):
#         ms2_fl = fl
#     if fl.endswith(".ms3"):
#         ms3_fl = fl



=======
>>>>>>> b57e85124cb1e518d151b77e429fa46303cb553f
# writeMs2ToMs3Info(ms2TriggerMS3dic, "ms2_ms3_linage.txt")

# def load_mass_table()
mH = 1.00782
mS = 31.9720
mH2O = 18.0105
mShort = 54.011
len_dic = {}


# 判断两个离子是否成对
def isPaired(a_list, b_list, tolppm=30):
    mz_a, z_a = a_list[:2]
    mz_b, z_b = b_list[:2]
    if z_a == z_b:
        if z_a == 0:
            return False, 50, 0
        else:
            min_mz = min(mz_a, mz_b)
            max_mz = max(mz_a, mz_b)
            ismatch, delta = compareTwoNum(min_mz+31.972/z_b, max_mz, tolppm)
            return ismatch, delta, z_b
    else:
        if min(z_a, z_b) > 0:
            return False, 50, 0
        else:
            z = max(z_a, z_b)
            min_mz = min(mz_a, mz_b)
            max_mz = max(mz_a, mz_b)
            ismatch, delta = compareTwoNum(min_mz+31.972/z, max_mz, tolppm)
            return ismatch, delta, z


# 去冗余，将字典里成对离子中的大质量的删掉，保留小的
def remove_large_in_pair(all_ms3):
    non_re_info = []
    i = 0
    while i < len(all_ms3)-1:
        for j in range(i+1, len(all_ms3)):
            ismatch, delta, z = isPaired(all_ms3[i], all_ms3[j])
            # print(ismatch, delta)
            if ismatch:
                if all_ms3[i] < all_ms3[j]:
                    a = list(all_ms3[i])
                else:
                    a = list(all_ms3[j])
                a[1] = z
                a.append(delta)
                non_re_info.append(tuple(a))
                del all_ms3[i]
                del all_ms3[j-1]
                break
        if not ismatch:
            non_re_info.append(all_ms3.pop(i))

    if all_ms3 != []:
        non_re_info.append(all_ms3[-1])
    return non_re_info


# 判断是否为alpha，beta的短的形式
def isSSFeatureIon(paired_ms3, ms3, ms2):
    mz1, z1 = paired_ms3[:2]
    mz2, z2 = ms3[:2]
    mpre, zpre = ms2[1:3]
    if z2 != 0:
        mzpre_from_ms3 = ((mz1- mH)*z1 + (mz2-mH)*z2 + mS + mH2O + mH * zpre)/zpre
        ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
        return ismatch, delta
    else:
        for z in range(1, zpre):
            mzpre_from_ms3 = ((mz1- mH)*z1 + (mz2-mH)*z2 + mS + mH2O + mH * zpre)/zpre
            ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
            if ismatch:
                return ismatch, delta
    
    return False, 50


# 判断是否为alpha短，beta长的形式
def isSLFeatureIon(paired_ms3, ms3, ms2):
    mz1, z1 = paired_ms3[:2]
    mz2, z2 = ms3[:2]
    mpre, zpre = ms2[1:3]
    if z2 != 0:
        mzpre_from_ms3 = ((mz1- mH)*z1 + (mz2-mH)*z2 + mH2O + mH * zpre)/zpre
        ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
        return ismatch, delta
    else:
        for z in range(1, zpre):
            mzpre_from_ms3 = ((mz1- mH)*z1 + (mz2-mH)*z + mH2O + mH * zpre)/zpre
            ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
            if ismatch:
                return ismatch, delta
        return False, 50
    

# 判断是否含有成对的信息
def isContainPair(non_re_info):
    paired_ms3 = [x for x in non_re_info if len(x) == 4]
    if len(paired_ms3) > 0:
        return True, paired_ms3[0]
    else:
        return False, None


# 判断所有ms3里是否有和确定的为alpha，beta特征离子
def judge_alpha_beta(all_ms3, paired_ms3, ms2):
    for ms3 in all_ms3:
        if ms3 != paired_ms3:
            isSS, delta_ss = isSSFeatureIon(paired_ms3, ms3, ms2)
            if isSS:
                return True, "isSS", [paired_ms3, ms3], delta_ss
            else:
                isSL, delta_sl = isSLFeatureIon(paired_ms3, ms3, ms2)
                if isSL:
                    return True, "isSL", [paired_ms3, ms3], delta_sl
    return False, None, None, 50


# 判断是否为短的mono形式
def isSMono_Feature(ms3, ms2):
    mz1, z1 = ms3[:2]
    mpre, zpre = ms2[1:3]
    if z1 != 0:
        mzpre_from_ms3 = ((mz1-mH) * z1 + mS + mH2O*2 + mShort+ mH * zpre)/zpre
        ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
        return ismatch, delta
    else:
        for z in range(1, zpre):
            mzpre_from_ms3 = ((mz1-mH) * z + mS + mH2O*2 + mShort+ mH * zpre)/zpre
            ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
            if ismatch:
                return ismatch, delta
        return False, 50


#判断是否为长的mono形式
def isLMono_Feature(ms3, ms2):
    mz1, z1 = ms3[:2]
    mpre, zpre = ms2[1:3]
    if z1 != 0:
        mzpre_from_ms3 = ((mz1-mH) * z1 + mH2O*2 + mShort+ mH * zpre)/zpre
        ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
        return ismatch, delta
    else:
        for z in range(1, zpre):
            mzpre_from_ms3 = ((mz1-mH) * z + mH2O*2 + mShort+ mH * zpre)/zpre
            ismatch, delta = compareTwoNum(mpre, mzpre_from_ms3, 20)
            if ismatch:
                return ismatch, delta
        return False, 50


# 判断是否含有mono
def judge_mono(non_re_info, ms2):
    if non_re_info == []:
        return False, None, None, 50
    else:
        for ms3 in non_re_info:
            isSM, delta_sm = isSMono_Feature(ms3, ms2)
            if isSM:
                return True, 'isSM', ms3, delta_sm
            else:
                isLM, delta_lm = isLMono_Feature(ms3, ms2)
                if isLM:
                    return True, 'isLM', ms3, delta_lm

        return False, None, None, 50


<<<<<<< HEAD
def main_flow(raw_name):
    ms2_fl = raw_name[:-4] + ".ms2"
    ms3_fl = raw_name[:-4] + ".ms3"
    if not os.path.exists(os.path.join(os.getcwd(), ms2_fl)) or \
        not os.path.exists(os.path.join(os.getcwd(), ms3_fl)):
        print("Check your files")
    else:
        ms3InfoDic = readms3Info(ms3_fl)
        numTotalMS3 = len(ms3InfoDic)
        ms2TriggerMS3dic = getLinageInfo(ms3InfoDic)
        totalms2, ms2TriggerMS3dic = readms2Info(ms2_fl, ms2TriggerMS3dic)
        xl = open(raw_name[:-4] + "_xlink.txt", 'w')
        mn = open(raw_name[:-4] + "_mono.txt", 'w')
        both_ms2_num = 0
        mono_num = 0
        for scan in ms2TriggerMS3dic:
            ms2_info = ms2TriggerMS3dic[scan]
            ms2 = ms2_info[0]
            all_ms3 = ms2_info[1:]
            print(ms2)
            print(all_ms3)
            non_re_info = remove_large_in_pair(all_ms3)
            print(non_re_info)
            contain_pair, paired_ms3 = isContainPair(non_re_info)
            if contain_pair:
                if len(non_re_info) >= 2:
                    isAlBe, isType, alpha_beta, delta_xl = judge_alpha_beta(non_re_info, paired_ms3, ms2)
                    if isAlBe:
                        print(scan)
                        print(isType)
                        print(alpha_beta, ms2)
                        wlist = [ms2, isType, delta_xl]
                        wlist.extend(alpha_beta)
                        xl.write("\t".join([str(x) for x in wlist])+"\n")
                        both_ms2_num += 1
                        non_re_left = [x for x in non_re_info if x not in alpha_beta]
                        isMono, mType, matched_ms3, delta_mn = judge_mono(non_re_left, ms2)
                        if isMono:
                            print(mType)
                            wlist = [ms2, matched_ms3, mType, delta_mn]
                            mn.write("\t".join([str(x) for x in wlist])+"\n")
                            mono_num += 1
                    else:
                        isMono, mType, matched_ms3, delta_mn = judge_mono(non_re_info, ms2)
                        if isMono:
                            print(scan)
                            print(mType)
                            wlist = [ms2, matched_ms3, mType, delta_mn]
                            mn.write("\t".join([str(x) for x in wlist])+"\n")
                            mono_num += 1
=======
def ms2_ms3_dsso_report(ms2_path, ms3_path):
    ms3InfoDic = readms3Info(ms3_path)
    numTotalMS3 = len(ms3InfoDic)
    ms2TriggerMS3dic = getLinageInfo(ms3InfoDic)
    totalms2, ms2TriggerMS3dic = readms2Info(ms2_path, ms2TriggerMS3dic)
    xl = open("xlink.txt", 'w')
    mn = open("mono.txt", 'w')
    both_ms2_num = 0
    mono_num = 0
    for scan in ms2TriggerMS3dic:
        ms2_info = ms2TriggerMS3dic[scan]
        ms2 = ms2_info[0]
        all_ms3 = ms2_info[1:]
        # print(ms2)
        # print(all_ms3)
        non_re_info = remove_large_in_pair(all_ms3)
        # print(non_re_info)
        contain_pair, paired_ms3 = isContainPair(non_re_info)
        if contain_pair:
            if len(non_re_info) >= 2:
                isAlBe, isType, alpha_beta, delta_xl = judge_alpha_beta(non_re_info, paired_ms3, ms2)
                if isAlBe:
                    print(scan)
                    print(isType)
                    print(alpha_beta, ms2)
                    wlist = [ms2, isType, delta_xl]
                    wlist.extend(alpha_beta)
                    xl.write("\t".join([str(x) for x in wlist])+"\n")
                    both_ms2_num += 1
                    non_re_left = [x for x in non_re_info if x not in alpha_beta]
                    isMono, mType, matched_ms3, delta_mn = judge_mono(non_re_left, ms2)
                    if isMono:
                        print(mType)
                        wlist = [ms2, matched_ms3, mType, delta_mn]
                        mn.write("\t".join([str(x) for x in wlist])+"\n")
                        mono_num += 1
>>>>>>> b57e85124cb1e518d151b77e429fa46303cb553f
                else:
                    isMono, mType, matched_ms3, delta_mn = judge_mono(non_re_info, ms2)
                    if isMono:
                        print(scan)
                        print(mType)
                        wlist = [ms2, matched_ms3, mType, delta_mn]
                        mn.write("\t".join([str(x) for x in wlist])+"\n")
                        mono_num += 1
<<<<<<< HEAD
        xl.close()
        mn.close()
        print(both_ms2_num, mono_num)
        print(both_ms2_num/len(ms2TriggerMS3dic))
        print(mono_num/len(ms2TriggerMS3dic))
        return both_ms2_num, mono_num, both_ms2_num/len(ms2TriggerMS3dic), mono_num/len(ms2TriggerMS3dic)

if __name__ == "__main__":
    b = open('report.csv', 'w')
    for fl in os.listdir(wk_dir):
        if fl.endswith(".raw"):
            both_ms2_num, mono_num, xlink_ratio, mono_ratio = main_flow(fl)
            b.write("%s,%d, %d, %f,%f\n" % (fl, both_ms2_num, mono_num, xlink_ratio, mono_ratio))
    b.close()


=======
            else:
                isMono, mType, matched_ms3, delta_mn = judge_mono(non_re_info, ms2)
                if isMono:
                    print(scan)
                    print(mType)
                    wlist = [ms2, matched_ms3, mType, delta_mn]
                    mn.write("\t".join([str(x) for x in wlist])+"\n")
                    mono_num += 1
    xl.close()
    mn.close()
    num_triggered_ms2 = len(ms2TriggerMS3dic)
    print(both_ms2_num, mono_num)
    print(both_ms2_num/num_triggered_ms2)
    print(mono_num/num_triggered_ms2)


if __name__ == "__main__":
    ms2_ms3_dsso_report('./DSSO_CID_MS2_CID_MS3_T1.ms2', './DSSO_CID_MS2_CID_MS3_T1.ms3')
>>>>>>> b57e85124cb1e518d151b77e429fa46303cb553f
