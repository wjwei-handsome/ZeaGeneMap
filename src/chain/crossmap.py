import sys
import pandas as pd
import pybedtools
from src.chain import chain
from tqdm import tqdm

from ..utils.utils import search_file
from src.utils.logger import logger

tqdm.pandas(desc='Processing', colour='green')

NULL_REGION = 'chrnull\t0\t0\t+'


def get_region(x, min_ratio: float = 0.8):
    chrom = x['chrom']
    start = x['start']
    end = x['end']
    strand = x['strand']
    # genename = x['genename']
    matches = chain.map_coordinates(bx_ivl, chrom, start, end, strand)
    if (matches is None) or (len(matches) % 2 != 0):
        return NULL_REGION
    # when matches == 2, there is one-to-one match (i.e. 100% match)
    if len(matches) == 2:
        region = '\t'.join([str(_) for _ in matches[1]])
        return region
    if len(matches) > 2:
        query_m = matches[::2]
        query_m_nt = sum([i[2] - i[1] for i in query_m])  # sum([3,2])
        target_m = matches[1::2]
        target_m_chroms = set([i[0] for i in target_m])
        target_m_starts = [i[1] for i in target_m]
        target_m_ends = [i[2] for i in target_m]
        target_m_strand = set([i[3] for i in target_m]).pop()
        map_ratio = query_m_nt / (end - start)
        if map_ratio >= min_ratio:
            if len(target_m_chroms) == 1:
                target_m_chrom = target_m_chroms.pop()
                target_m_start = min(target_m_starts)
                target_m_end = max(target_m_ends)
                region = '\t'.join([target_m_chrom, str(
                    target_m_start), str(target_m_end), target_m_strand])
                return region
            else:
                return NULL_REGION
        else:
            return NULL_REGION


# def string_bed_ivl(x, min_ovp: float = 0.8):
#     if x != 'None':
#         str_bed = pybedtools.bedtool.BedTool(x, from_string=True)
#         intersect_ivl = target_bed.intersect(
#             str_bed, nonamecheck=True, F=min_ovp, sorted=True)
#         if intersect_ivl:
#             result = []
#             intersect_ivl_lines = intersect_ivl.__str__().splitlines()
#             for line in intersect_ivl_lines:
#                 line_chr, line_start, line_end, line_junk, line_name, line_strand = line.split(
#                     '\t')
#                 result.append(
#                     f'{line_name}@({line_chr}:{line_start}-{line_end}-{line_strand})')
#             str_result = ','.join(result)
#             return str_result
#         else:
#             result = []
#             for ivl in str_bed:
#                 result.append(
#                     f'({ivl.chrom}:{ivl.start}-{ivl.end})')
#             str_result = ','.join(result)
#             return str_result
#     else:
#         return 'None'

def get_intersect(region_df: pd.DataFrame, overlap: float = 0.8) -> pybedtools.bedtool.BedTool:
    """
    Get intersect between target and query bed files.
    """
    region_list = region_df.to_list()
    region_bed_str = '\n'.join(region_list)
    region_bed = pybedtools.bedtool.BedTool(region_bed_str, from_string=True)
    intersect_bed = region_bed.intersect(target_bed, nonamecheck=True,
                                         f=overlap, wo=True, loj=True)
    return intersect_bed


def get_intersect_df(intersect_bed: pybedtools.bedtool.BedTool) -> pd.DataFrame:
    """_summary: Get intersect df from a bed.

    Args:
        intersect_bed (pybedtools.bedtool.BedTool): _description_

    Returns:
        pd.DataFrame: _description_
    """
    lst = []
    for ivl in intersect_bed:
        chrom = ivl.chrom
        region_start, region_end = ivl.start, ivl.end
        gene_start, gene_end = ivl.fields[5:7]
        if chrom == 'chrnull':
            lst.append('NoneRegion')
        else:
            if gene_start == '-1' or gene_end == '-1':
                lst.append(f'{chrom}:{region_start}-{region_end}')
            else:
                target_gene = ivl.fields[7]
                lst.append(
                    f'{target_gene}@{chrom}:{max([region_start,int(gene_start)])}-{min([region_end,int(gene_end)])}')
    df = pd.DataFrame(lst, columns=['crossmap'])
    return df


def get_crossmap_df(raw_df: pd.DataFrame,
                    query_g: str,
                    target_g: str,
                    work_dir: str,
                    ratio: float = 0.8,
                    overlap: float = 0.8):
    chain_fp = search_file(work_dir, 'chain.gz',
                           query_g=query_g, target_g=target_g, type='single')
    logger.info('Start to combine crossmap evidence')
    logger.info(f"Searching Files:\n chain file: {chain_fp}")
    logger.info(f"Reading chain file: {chain_fp}")
    global bx_ivl
    bx_ivl = chain.read_chain_file(chain_fp)[0]
    q_bed_fp = search_file(work_dir, 'bed', single_g=query_g, type='single')
    logger.info(f"Searching Files:\n query genome bed file:{q_bed_fp}")
    t_bed_fp = search_file(work_dir, 'bed', single_g=target_g, type='single')
    logger.info(f"Searching Files:\n taget genome bed file: {t_bed_fp}")
    q_bed_df = pd.read_csv(q_bed_fp, sep='\t', header=None)
    q_bed_df.columns = ['chrom', 'start', 'end', 'genename', 'junk', 'strand']
    raw_bed_df = pd.merge(raw_df, q_bed_df, on='genename',
                          how='left').drop(columns=['junk'])
    global target_bed
    target_bed = pybedtools.BedTool(t_bed_fp).sort()
    raw_bed_df['region'] = raw_bed_df.progress_apply(
        get_region, axis=1, args=(ratio,))

    intersect_bed = get_intersect(raw_bed_df['region'], overlap)
    intersect_df = get_intersect_df(intersect_bed)

    # concat method
    if len(intersect_df) == len(raw_bed_df):
        return pd.concat([raw_bed_df['genename'], intersect_df], axis=1)
    else:
        logger.error('intersect_df and raw_bed_df have different length')
        sys.exit(1)

    # raw_bed_df['crossmap'] = raw_bed_df['region'].progress_apply(
    #     string_bed_ivl, args=(overlap,))
    # return raw_bed_df[['genename', 'crossmap']]
