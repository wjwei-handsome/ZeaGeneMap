import pandas as pd
from src.utils.logger import logger
from ..utils.utils import search_file


# def get_rbh_result(gene_name, rbh_df):
#     result = []
#     filter_df = rbh_df.query(f'qseqid_x=="{gene_name}"')
#     if filter_df.empty:
#         return 'None'
#     else:
#         for index, row in filter_df.iterrows():
#             result.append(row['qseqid_y'])
#         result_str = ','.join(result)
#         return result_str


def get_rbh_df(query_g: str, target_g: str, work_dir: str, e_filter: float = 1e-10) -> pd.DataFrame:
    q_t_blast_fp, t_q_blast_fp = search_file(
        work_dir, 'blast', query_g, target_g, type='rec')
    logger.info('Start to combine rbh evidence')
    logger.info(
        f'Searching Files:\n query -> target blast file: {q_t_blast_fp}\n target -> query blast file: {t_q_blast_fp}')
    fwd_df = pd.read_csv(q_t_blast_fp, sep='\t', header=None)
    rev_df = pd.read_csv(t_q_blast_fp, sep='\t', header=None)
    column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                    'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    fwd_df.columns = column_names
    rev_df.columns = column_names
    fwd_df = fwd_df.drop_duplicates(subset=['qseqid'], keep='first')
    rev_df = rev_df.drop_duplicates(subset=['qseqid'], keep='first')
    rbh = pd.merge(fwd_df,
                   rev_df[['qseqid', 'sseqid']],
                   left_on='sseqid',
                   right_on='qseqid',
                   how='outer')
    rbh = rbh.loc[rbh.qseqid_x == rbh.sseqid_y]
    rbh = rbh.groupby(['qseqid_x', 'sseqid_x']).max().reset_index()[
        ['qseqid_x', 'sseqid_x']]
    rbh.columns = ['genename', 'rbh']
    return rbh
