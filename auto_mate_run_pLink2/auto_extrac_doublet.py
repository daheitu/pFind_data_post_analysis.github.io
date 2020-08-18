import os
import time
os.chdir(r"C:\Users\Yong Cao\Documents\github\pLink_DSSO_output_doublet\pLink_DSSO_output_doublet\pLink\bin")
from automate_pLink2 import count_pf2_file, makedir

def run_searcher(raw_path, spec_type, mgf_name):
    f = open("./plink2_plink_dsso.plink").readlines()
    output_name = "%s_output" % mgf_name[:-4]
    makedir(raw_path, output_name)
    plink_para_name = "plink2.plink" 
    plink_path = os.path.join(raw_path, output_name, plink_para_name)
    b = open(plink_path, 'w')
    # para_massF, flow_type_inner = type_para_dic[flow_type]
    for line in f:
        if line.startswith("spec_num"):
            b.write(line)
            break
        else:
            # if line.startswith("processor_num ="):
            #     b.write("processor_num = 3\n")
            if line.startswith("result_output_path"):
                b.write("result_output_path = " + os.path.join(raw_path, output_name) +"\n")
            # elif line.startswith("db_name"):
            #     b.write("db_name = " + db_name + "\n")
            # elif line.startswith("evalue ="):
            #     b.write("evalue = 1\n")
            # elif line.startswith("flow_type"):
            #     b.write("flow_type = %s\n" % flow_type_inner)
            # elif line.startswith("mass_filter_4_ion_index_flow"):
            #     b.write("mass_filter_4_ion_index_flow = %s\n" % para_massF)
            # elif line.startswith("linker1"):
            #     b.write("linker1 = %s\n" % linker)
            else:
                b.write(line)
        
    cur_time = time.strftime("%Y-%m-%d ", time.gmtime())
    b.write("spec_title = " + "_"+ cur_time+"\n")
    b.write("spec_type = %s\n" % spec_type)
    # fl_list = count_pf2_file(raw_path, spec_type)
    # for line in fl_list:
    b.write("spec_path1 = %s" % (os.path.join(raw_path, mgf_name))+"\n")
    b.write("[quant]\n")
    b.write("quant = 1|None\n")
    b.close()
    cmd_search = "searcher.exe " + plink_path
    os.system(cmd_search)


if __name__ == "__main__":
    raw_path = r"F:\data_from_paper\maxlinker_MCP\JinLiang10366422\raw"
    spec_type = "mgf"
    for fl in os.listdir(raw_path):
        if fl.endswith(spec_type):
            run_searcher(raw_path, spec_type, fl)