# coding = utf-8

import os, time

raw_p = r"G:\msData\20200419\BSA"
plink_bin_path = r"E:\pFindStudio\pLink\2.3.9_20200327\bin"
plink_para_demo = r"E:\pFindStudio\pLink_test_20200330\pLink_test_20200330\2.pLink2_plus_score_massfilter\2.pLink2_plus_score_massfilter.plink"
db_name = "bsa_con" # "CNGP_con" # "GST_sequence_con" # "BSA_with_resoure"

#flow_type = 0 # 0 for +score； 1 for +mass filter； 2 for pLink-DSSO


type_name_dic = {0:"_score", 1:"_mass_filter", 2:"_plink_dsso"}

type_para_dic = {0:('0', "FT_ION_INDEX"),
                1:('1', "FT_ION_INDEX"),
                2:('1', 'FT_SteppedCleavage_PEP_INDEX')}


def makedir(rootpath, dir_name):
    path = os.path.join(rootpath, dir_name)
    if os.path.exists(path):
        pass
    else:
        os.makedirs(path)


def get_pParse_para(raw_path):
    f = open(os.path.join(plink_bin_path, "pParse.para")).readlines()
    makedir(raw_path, "pParse_para")
    pParse_path = os.path.join(raw_path, "pParse_para", "pParse.para")
    b = open(pParse_path, 'w')
    for line in f:
        if line.startswith("datapath ="):
            b.write("datapath = " + raw_path +"\n")
        else:
            b.write(line)
    b.close()
    os.chdir(plink_bin_path)
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


def run_searcher(flow_type, linker, raw_path):
    f = open(plink_para_demo).readlines()
    output_name = "output%s" % (type_name_dic[flow_type])
    makedir(raw_path, output_name)
    plink_para_name = "plink2%s.plink" % (type_name_dic[flow_type])
    plink_path = os.path.join(raw_path, output_name, plink_para_name)
    b = open(plink_path, 'w')
    para_massF, flow_type_inner = type_para_dic[flow_type]
    for line in f:
        if line.startswith("[spectrum]"):
            b.write(line)
            break
        else:
            if line.startswith("processor_num ="):
                b.write("processor_num = 3\n")
            elif line.startswith("result_output_path"):
                b.write("result_output_path = " + os.path.join(raw_path, output_name) +"\n")
            elif line.startswith("db_name"):
                b.write("db_name = " + db_name + "\n")
            elif line.startswith("evalue ="):
                b.write("evalue = 1\n")
            elif line.startswith("flow_type"):
                b.write("flow_type = %s\n" % flow_type_inner)
            elif line.startswith("mass_filter_4_ion_index_flow"):
                b.write("mass_filter_4_ion_index_flow = %s\n" % para_massF)
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


def main_flow(raw_path, linker):
    get_pParse_para(raw_path)
    for i in range(3):    
        flow_type = i
        run_searcher(flow_type, linker, raw_path)        


if __name__ == "__main__":
    for root, dirs, fls in os.walk(raw_p):
        if root.split("\\")[-1] == "DSSO":
            print(root)
            raw_path = root
            linker = "DSSO"
            main_flow(raw_path, linker)