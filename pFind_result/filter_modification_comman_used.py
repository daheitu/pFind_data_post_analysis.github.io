import os
# from typing import Sized

wkdir = r"Z:\pFind_work_space\Juri_DSS\result" #important! the result path of pFind tastk 
pro_looking = ""#"wins_seg"  # important! enter the protein of interset, if no name filled, the program will serach all identified proteins  
modify = "Xlink_DSS[156][K]"#"Carbamidomethyl[C]"  # important! enter the modification type you want search 


##############################Don't change the following lines###########################
def find_spectra(pep,mod, num, spec):
    spec_title_list = []
    score_list = []
    for line in spec:
        linelist = line.split("\t")
        if len(linelist) == 19:
            ps = linelist[5]
            ms = linelist[10]
            if ps == pep and ms == mod:
                spec_title_list.append(linelist[0])
                score_list.append(float(linelist[9]))
    return spec_title_list[:num], score_list[:num]


#获取raw文件的列表
def get_raw_file_list():
    raw_list = []
    spec_f = open(r"../param/pFind.cfg").readlines()
    for line in spec_f:
        if line.startswith("msms") and line.endswith(".pf2\n"):
            raw = line.split("\\")[-1][:-11]
            raw_list.append(raw)
    print(raw_list)
    return raw_list


def get_all_id_proteins():
    proteins = []
    f = open("pFind.protein").readlines()
    for line in f:
        linelist = line.split("\t")
        if linelist[0].isdigit():
            protein = linelist[1]
            if protein not in proteins:
                proteins.append(protein)
        elif linelist[0] == "" and linelist[1] in ["SubSet", "SameSet"]:
            protein = linelist[2]
            if protein not in proteins:
                proteins.append(protein)
    print(proteins)
    return proteins


def get_modi_info(pro_looking, modify):
    modi_dic = {}
    f = open("pFind.protein").readlines()
    spec_f = open("pFind.spectra").readlines()

    for line in f[:]:
        linelist = line.split("\t")
        if len(linelist) == 19:        
            pros = linelist[10]
            mods = linelist[8]
            if pro_looking in pros and modify in mods:
                pep = linelist[3]
                pros_list = pros[:-1].split('/')
                pep_pos_list = linelist[11][:-1].split('/')
                mod_list = mods[:-1].split(";")
                spec_num = int(linelist[-1])
                score = float(linelist[7])
                spec = linelist[-3]
                for mod_info in mod_list:
                    md_site, md = mod_info.split(',')
                    if md  == modify:
                        for i in range(len(pros_list)):
                            pro = pros_list[i]
                            if pro == pro_looking:
                                pep_pos = int(pep_pos_list[i].split(',')[0])
                                mod_pos = pep_pos + int(md_site)
                                titles, scores = find_spectra(pep, mods, spec_num, spec_f)
                                if mod_pos not in modi_dic:
                                    modi_dic[mod_pos] = [spec_num, scores,[pep], titles]
                                else:
                                    modi_dic[mod_pos][0] += spec_num
                                    modi_dic[mod_pos][1].extend(scores)
                                    modi_dic[mod_pos][2].append(pep)
                                    # modi_dic[mod_pos][3].append(spec)
                                    modi_dic[mod_pos][3].extend(titles)
    
    # print(modi_dic)
    return modi_dic

def find_top3_spec(scores, specrums):
    if len(scores) != len(specrums):
        print("wrong")
    else:
        info_list = []
        for i in range(len(scores)):
            info_list.append((scores[i], specrums[i]))
        info_list.sort(key=lambda x:x[0], reverse= True)
        return [x[1] for x in info_list][:3]


def format_dic(modi_dic):
    new_dic = {}
    spectrum_dic = {}
    for site in modi_dic:
        info_list = modi_dic[site]
        specs = info_list[-1]
        scores = info_list[-3]
        top3_spectrum = find_top3_spec(scores, specs)
        spectrum_dic[site] = top3_spectrum
        new_dic[site] = {}
        for i in range(len(specs)):
            raw = specs[i].split(".")[0]
            # print(raw)
            if raw not in new_dic[site]:
                new_dic[site][raw] = [1, scores[i]]
                
            else:
                new_dic[site][raw][0] += 1
                if scores[i] < new_dic[site][raw][1]:
                    new_dic[site][raw][1] = scores[i]
    # print(new_dic, spectrum_dic)
    return new_dic, spectrum_dic



def write_info_2_file(pro_looking, b, raw_list, modify):
    print(pro_looking)
    modi_dic = get_modi_info(pro_looking, modify)
    if modi_dic != {}:
        reformat_dic, specrum_dic = format_dic(modi_dic)
        for mod_pos in reformat_dic:
            wlist = [pro_looking, mod_pos, modify[-2], modify, modi_dic[mod_pos][0], min(modi_dic[mod_pos][1])]
            info_dic = reformat_dic[mod_pos]
            for raw in raw_list:
                if raw in info_dic:
                    wlist.extend(info_dic[raw])
                else:
                    wlist.extend(["", ""])
            # wlist = [mod_pos,"S", info_list[0], info_list[1], len(info_list[2]), info_list[3][0], info_list[4]]
            # b.write(pro_looking+"\n")
            wlist.extend(specrum_dic[mod_pos])
            print(wlist)
            b.write(",".join([str(ele) for ele in wlist])+"\n")



def main(wkdir, pro_looking, modify):
    os.chdir(wkdir)
    report_file_name = modify + "_"+ "_report.csv"
    b = open(report_file_name, 'w')
    titleline = ["Protein", "site","aa", "Modification", "total_spec", "best_score"]
    
    raw_list = get_raw_file_list()
    
    for raw in raw_list:
        titleline.append(raw+"_spec")
        titleline.append(raw+"_score")
    b.write(",".join(titleline)+"\n")

    if pro_looking !=  "":
        write_info_2_file(pro_looking, b, raw_list, modify)
    else:
        proteins = get_all_id_proteins()
        for prot in proteins:
            write_info_2_file(prot, b , raw_list, modify)

    b.close()


if __name__ == "__main__":
    main(wkdir, pro_looking, modify)
    print("Well Done!")