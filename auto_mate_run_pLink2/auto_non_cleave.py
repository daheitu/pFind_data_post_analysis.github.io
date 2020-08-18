# coding = utf-8
# author: Yong Cao

# coding = utf-8

import os, time
from automate_pLink2 import is_searched


raw_p = r"F:\MS_DATA_STORAGE\pku_RH\DSS"
plink_bin_path = r"E:\pFindStudio\pLink2.3.9_0415\bin"
# plink_para_demo = r"E:\pFindStudio\pLink_test_20200330\pLink_test_20200330\2.pLink2_plus_score_massfilter\2.pLink2_plus_score_massfilter.plink"
# db_name = "synthetic_pep_con" #"Lactoferrin_con" # "CNGP_con" # "GST_sequence_con" # "BSA_with_resoure"

db_name_dic = {
    "hGBP":"lipeng_con",
    "Lysozyme":"lysozyme_con",
    "PUD12":"PUD12_con",
    "UTPA":"UTPA_sub complex_con"
}
#                 "BSA":"bsa_con", 
#                 "CNGP": "CNGP_con",
#                 "GST": "GST_sequence_con",
#                 "Lacto": "Lactoferrin_con"}


def makedir(rootpath, dir_name):
    path = os.path.join(rootpath, dir_name)
    if os.path.exists(path):
        pass
    else:
        os.makedirs(path)


def get_pParse_para(raw_path):
    os.chdir(plink_bin_path)
    f = open("pParse.para").readlines()
    makedir(raw_path, "pParse_para")
    pParse_path = os.path.join(raw_path, "pParse_para", "pParse.para")
    b = open(pParse_path, 'w')
    for line in f:
        if line.startswith("datapath ="):
            b.write("datapath = " + raw_path +"\n")
        else:
            b.write(line)
    b.close()
    cmd_pparse = "pParse.exe " + pParse_path
    os.system(cmd_pparse)


def count_pf2_file(raw_path):
    fl_list = []
    i  = 0
    for fl in os.listdir(raw_path):
        if fl.endswith("pf2"):
            print(fl)
            i += 1
            fl_list.append("spec_path%d = %s" % (i, os.path.join(raw_path, fl)))
    fl_list.insert(0, "spec_num = %d" % i)
    return fl_list

#print(count_pf2_file(raw_path))


def run_searcher(linker, raw_path, db_name):
    f = open("./demo_DSS.plink").readlines()
    output_name = "output"
    makedir(raw_path, output_name)
    plink_para_name = "plink2.plink" 
    plink_path = os.path.join(raw_path, output_name, plink_para_name)
    b = open(plink_path, 'w')
    for line in f:
        if line.startswith("[spectrum]"):
            b.write(line)
            break
        else:
            if line.startswith("result_output_path"):
                b.write("result_output_path = " + os.path.join(raw_path, output_name) +"\n")
            elif line.startswith("db_name"):
                b.write("db_name = " + db_name + "\n")
            elif line.startswith("evalue ="):
                b.write("evalue = 1\n")
            elif line.startswith("linker1"):
                b.write("linker1 = %s\n" % linker)
            else:
                b.write(line)
        
    cur_time = time.strftime("%Y-%m-%d ", time.gmtime())
    b.write("spec_title = " + db_name + "_"+ cur_time+"\n")
    b.write("spec_type = pf\n")
    fl_list = count_pf2_file(raw_path)
    for line in fl_list:
        b.write(line+"\n")
    b.write("[quant]\n")
    b.write("quant = 1|None\n")
    b.close()
    cmd_search = "searcher.exe " + plink_path
    os.system(cmd_search)


def main_flow(raw_path, linker, db_name):
    get_pParse_para(raw_path)
    run_searcher(linker, raw_path, db_name)


if __name__ == "__main__":
    main_flow(raw_p, "DSS", "X.laevis_NPL_200")
    # for root, dirs, fls in os.walk(raw_p):
    #     linker = root.split("\\")[-1]
    #     if linker in ["DSS", "BSMEG"]:
    #         pro_name = root.split("\\")[-2]
    #         raw_path = root
    #         # linker = "DSS"
    #         if not is_searched(raw_path):
    #             print(root)
    #             db_name = db_name_dic[pro_name]
    #             main_flow(raw_path, linker, db_name)
    #             print("sleep time: 10 mins")
    #             time.sleep(120)