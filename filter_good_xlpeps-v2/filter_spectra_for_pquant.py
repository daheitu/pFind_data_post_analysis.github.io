import os
# from generate_combination import load_peptide_dic
# from judge_xlpep_num import get_all_pep_dic, judge_peptide
# all_pep_dic = get_all_pep_dic()


# xl_spec_path = r"G:\msData\synthetic_pepteide_rawdata\CV3_CV7\DSSO\output_score\reports\synthetic_pep_con_2020-04-30.filtered_cross-linked_spectra.csv"

def split_xl_spec(xl_spec_path):
    print(os.path.basename(xl_spec_path), os.path.dirname(xl_spec_path))
    f= open(xl_spec_path).readlines()

    raw_lines_dic = {}
    for line in f[1:]:
        linelist = line.split(',')
        # pep = linelist[4]
        title = linelist[1]
        raw_name = title.split('.')[0]
        if raw_name not in raw_lines_dic:
            raw_lines_dic[raw_name] = [line]
        else:
            raw_lines_dic[raw_name].append(line)

    # print(len(raw_lines_dic))
    dir_name = os.path.dirname(xl_spec_path)
    if len(raw_lines_dic) > 1:
        for raw_name in raw_lines_dic:
            file_name = raw_name + "_cross-linked_spectra_filter.csv"
            rep_path = os.path.join(dir_name, file_name)
            b = open(rep_path, 'w')
            b.write(f[0])
            for line in raw_lines_dic[raw_name]:
                b.write(line)
            b.close()
    