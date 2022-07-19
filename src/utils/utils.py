import os
import glob


def search_file(work_dir: str,
                surfix: str,
                query_g: str = '', target_g: str = '',
                type: str = 'rec',
                single_g: str = ''):
    abs_dir_path = os.path.abspath(work_dir)
    if type == 'rec':
        forward_fp = f"{query_g}.{target_g}.{surfix}"
        reverse_fp = f"{target_g}.{query_g}.{surfix}"
        fw_search = glob.glob(f"{abs_dir_path}/**/{forward_fp}")
        rv_search = glob.glob(f"{abs_dir_path}/**/{reverse_fp}")
        if len(fw_search) == 1 and len(rv_search) == 1:
            return fw_search[0], rv_search[0]
        else:
            exit(1)
    elif type == 'single':
        single_fp = f"{single_g}.{surfix}"
        single_search = glob.glob(f"{abs_dir_path}/**/{single_fp}")
        if len(single_search) == 1:
            return single_search[0]
        else:
            exit(1)
