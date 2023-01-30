from hmac import digest
import os

wk_dir = r"Z:\pFind_work_spacejuri_mito_dss_open\result"
# digest_fasta_path = r"Z:\STY_PROJ\pFind_search_dataset\LiuFAN_DSS_HELA_FAIMS\UP000005640_Homo_sapiens_20190921_reference_con_filter.fasta_enzyme"
#"D:\FastaDatabase\UP000005640_Homo_sapiens_20190921_reference_con.fasta_enzyme"
# "Z:\STY_PROJ\pFind_search_dataset\LiuFAN_DSS_HELA_FAIMS\UP000005640_Homo_sapiens_20190921_reference_con_filter.fasta_enzyme"
# "D:\FastaDatabase\uniprot-proteome_ Escherichia coli (strain K12)_20210123_xuzhancong.fasta_enzyme"
#"

# tgt_mod_list = ["Leiker_CLV[334][K]", "Leiker_CLV[334][S]", "Leiker_CLV[334][T]", "Leiker_CLV[334][Y]",\
#      "Leiker_CLV[334][G]", "Leiker_CLV[334][V]", "Leiker_CLV[334][L]"]

tgt_mod_list = ["Xlink_DSS[156][K]", "Xlink_DSS[156][S]", "Xlink_DSS[156][T]", "Xlink_DSS[156][Y]", "Xlink_DSS[156][G]", "Xlink_DSS[156][V]", "Xlink_DSS[156][L]"]

rep_file_name = wk_dir.split("\\")[-2] + ".csv"

def stac_spec_mods():
    f = open(os.path.join(wk_dir, "pFind.protein")).readlines()

    def is_wanted(proteins):
        protein_list = proteins.split("/")[:-1]
        for pro in protein_list:
            if not pro.startswith("REV_"):
                return True
        return False

    # print(is_wanted("sp|P78009|tsf/"))

    rep_dic = {}
    for mod in tgt_mod_list:
        rep_dic[mod] = [0, 0]

    pep_set = set()

    for line in f[2:]:
        linelist = line.split("\t")

        if linelist[0] == "" and linelist[1] == "":
            # print(line)
            seq = linelist[3]
            mods = linelist[8]
            proteins = linelist[10]
            spec_num = int(linelist[-1].strip())
            # print(seq, mods, proteins, spec_num)
            if is_wanted(proteins) and  (seq, mods) not in pep_set:
                pep_set.add((seq, mods))
                mods_list = mods.split(";")[:-1]
                # print(mods_list)
                for mod in tgt_mod_list:
                    for mod_site in mods_list:
                        if mod in mod_site:
                            rep_dic[mod][0] += spec_num
                            rep_dic[mod][1] += 1

    return rep_dic


def get_condidate_site(fasta_path_digst, amino_acid):
    f = open(fasta_path_digst).readlines()
    pep_set = set()
    num = 0
    for line in f:
        if not line.startswith(">"):
            pep_set.add(line.strip())
    for pep in pep_set:
        num += pep[:-1].count(amino_acid)
    
    return num


def main():
    b = open(rep_file_name, 'w')
    b.write("mod, spec_num, cond_pep, spec_num/cond_pep, pep_num, pep_num/cond_pep\n")
    rep_dic = stac_spec_mods()
    for mod in rep_dic:
        print(mod)
        spec_num = rep_dic[mod][0]
        pep_num = rep_dic[mod][1]
        # cond_pep = get_condidate_site(digest_fasta_path, mod[-2])
        # b.write("%s,%d,%d,%f,%d,%f\n" % (mod, spec_num, cond_pep, spec_num/cond_pep, pep_num, pep_num/cond_pep))
        # cond_pep = get_condidate_site(digest_fasta_path, mod[-2])
        b.write("%s,%d,%d\n" % (mod, spec_num, pep_num))
    
    b.close()


main()
# print(get_condidate_site(digest_fasta_path, 'K'))