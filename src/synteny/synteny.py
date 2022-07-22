import pandas as pd

from ..utils.utils import search_file
from src.utils.logger import logger


def get_rec_syn_df(query_g: str, target_g: str, work_dir: str) -> pd.DataFrame:
    q_t_synteny_fp, t_q_synteny_fp = search_file(
        work_dir, 'anchor', query_g, target_g, type='rec')
    logger.info('Start to combine syteny evidence')
    logger.info(
        f'Searching Files:\n query -> target synteny file: {q_t_synteny_fp}\n target -> query synteny file: {t_q_synteny_fp}')
    exclude_fwd = [i for i, line in enumerate(
        open(q_t_synteny_fp)) if line.startswith('#')]
    exclude_rev = [i for i, line in enumerate(
        open(t_q_synteny_fp)) if line.startswith('#')]

    fwd_df = pd.read_csv(q_t_synteny_fp, sep='\t', header=None,
                         usecols=[0, 1], skiprows=exclude_fwd)
    rev_df = pd.read_csv(t_q_synteny_fp, sep='\t', header=None,
                         usecols=[0, 1], skiprows=exclude_rev)
    fwd_df.columns = ['query_g', 'target_g']
    rev_df.columns = ['target_g', 'query_g']

    rev_syn_df = pd.merge(fwd_df,
                          rev_df,
                          on='query_g',
                          how='inner')
    rev_syn_df = rev_syn_df[rev_syn_df['target_g_x'] ==
                            rev_syn_df['target_g_y']].reset_index(drop=True)
    rev_syn_df.drop(columns=['target_g_y'], inplace=True)
    rev_syn_df.columns = ['genename', 'synteny']
    rev_syn_df = rev_syn_df.groupby('genename', as_index=False).agg(
        lambda x: ','.join(list(x)) if len(x) > 1 else x)
    return rev_syn_df
