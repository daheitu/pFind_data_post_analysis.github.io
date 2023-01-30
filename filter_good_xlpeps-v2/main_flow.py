import os
from pquant import gene_run_pquant
from cal_RT_time import cal_write_result
from filter_spectra_for_pquant import split_xl_spec

wk_dir = r'M:\TanDan_8std_dilute\BS3_INPUT\10'
plink_bin_path = r"E:\pFindStudio\pLink2.3.9_0415\bin"


def split_spec(plink_id_folder):
    for root, drs, fls in os.walk(os.path.join(plink_id_folder, "reports")):
        for fl in fls:
            if fl.endswith("filtered_cross-linked_spectra.csv"):
                split_xl_spec(os.path.join(root, fl))
        


def main_flow(raw_path):
    plink_id_folder = os.path.join(raw_path, "PrMix_output")
    gene_run_pquant(plink_bin_path, raw_path)
    
    split_spec(plink_id_folder)
    cal_write_result(plink_id_folder)


if __name__ == "__main__":
    # for root, drs, fls in os.walk(wk_dir):
    #     if root.split('\\')[-1] == "DSSO" and "CV3_CV7" not in root:
    #         print(root)
    #         main_flow(root)
    raw_path = r"M:\TanDan_8std_dilute\BS3_INPUT\1100"
    main_flow(raw_path)