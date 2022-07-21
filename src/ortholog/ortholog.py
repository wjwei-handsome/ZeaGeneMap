import pandas as pd
from ..utils.utils import search_file


def map_str_list(x: str) -> list:
    lst = x.split(',')
    lst = [i.strip() for i in lst]
    lst = [i.split(' ')[1] if len(i.split(' ')) == 2 else i for i in lst]
    if len(lst) == 1:
        return lst[0]
    return lst


def get_ortholog_df(query_g: str, target_g: str, work_dir: str):
    """
    get ortholog df
    """
    ortholog_fp = search_file(work_dir, 'tsv', type='ortholog')
    print(f"find:\n ortholog_fp: {ortholog_fp}")
    used_cols = [query_g, target_g]
    ortholog_df = pd.read_csv(ortholog_fp, sep='\t', header=0,
                              usecols=used_cols)[used_cols]
    ortholog_df.fillna('None', inplace=True)
    ortholog_df[query_g] = ortholog_df[query_g].map(map_str_list)
    ortholog_df[target_g] = ortholog_df[target_g].map(map_str_list)
    ortholog_df = ortholog_df.explode(query_g)
    ortholog_df = ortholog_df.reset_index(drop=True)
    ortholog_df.columns = ['genename', 'ortholog']
    ortholog_df['ortholog'] = ortholog_df['ortholog'].map(
        lambda x: ','.join(x) if isinstance(x, list) else x)
    # ortholog_df = ortholog_df.groupby('genename', as_index=False).agg(
    #     lambda x: list(x) if len(x) > 1 else x)

    return ortholog_df
