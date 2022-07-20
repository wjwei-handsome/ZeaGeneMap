import pandas as pd
import pybedtools
from src.chain import chain
from tqdm import tqdm

from ..utils.utils import search_file

tqdm.pandas()

# tmp_fp = 'tests/Zm-B73v5.Zv-TEO01.chain.gz'

# q_bed_fp = 'tests/Zm-B73v5.bed'
# t_bed_fq = 'tests/Zv-TEO01.bed'

# # raw_df = pd.read_csv('tests/Zm-B73v5.genelist', sep='\t', header=None)
# # raw_df.columns = ['genename']

# q_bed_df = pd.read_csv(q_bed_fp, sep='\t', header=None)
# q_bed_df.columns = ['chrom', 'start', 'end', 'genename', 'junk', 'strand']

# raw_bed_df = pd.merge(raw_df, q_bed_df, on='genename',
#                       how='left').drop(columns=['junk'])

# bx_ivl = chain.read_chain_file(tmp_fp)[0]


# target_bed = pybedtools.BedTool(t_bed_fq)


def get_region(x, min_ratio: float = 0.85):
    chrom = x['chrom']
    start = x['start']
    end = x['end']
    strand = x['strand']
    # genename = x['genename']
    matches = chain.map_coordinates(bx_ivl, chrom, start, end, strand)
    if (matches is None) or (len(matches) % 2 != 0):
        return 'None'
    # when matches == 2, there is one-to-one match (i.e. 100% match)
    if len(matches) == 2:
        region = '\t'.join([str(_) for _ in matches[1]])
        return region
    if len(matches) > 2:
        query_m = matches[::2]
        query_m_nt = sum([i[2]-i[1] for i in query_m])  # sum([3,2])
        # ODDS: [('chr1', 248908207, 248908210, '+'), ('chr1', 249058210, 249058212, '+')]
        target_m = matches[1::2]
        target_m_chroms = set([i[0] for i in target_m])
        target_m_starts = [i[1] for i in target_m]
        target_m_ends = [i[2] for i in target_m]
        target_m_strand = set([i[3] for i in target_m]).pop()
        #print (a_target_ends)
        map_ratio = query_m_nt/(end-start)
        if map_ratio >= min_ratio:
            if len(target_m_chroms) == 1:
                target_m_chrom = target_m_chroms.pop()
                target_m_start = min(target_m_starts)
                target_m_end = max(target_m_ends)
                region = '\t'.join([target_m_chrom, str(
                    target_m_start), str(target_m_end), target_m_strand])
                return region
            else:
                return 'None'
        else:
            return 'None'


def string_bed_ivl(x, min_ovp: float = 0.8):
    if x != 'None':
        str_bed = pybedtools.bedtool.BedTool(x, from_string=True)
        intersect_ivl = target_bed.intersect(
            str_bed, nonamecheck=True, F=min_ovp)
        if intersect_ivl:
            result = []
            intersect_ivl_lines = intersect_ivl.__str__().splitlines()
            for line in intersect_ivl_lines:
                line_chr, line_start, line_end, line_junk, line_name, line_strand = line.split(
                    '\t')
                result.append(
                    f'{line_name}@({line_chr}:{line_start}-{line_end}-{line_strand})')
            str_result = ','.join(result)
            return str_result
        else:
            result = []
            for ivl in str_bed:
                result.append(
                    f'({ivl.chrom}:{ivl.start}-{ivl.end})')
            str_result = ','.join(result)
            return str_result
    else:
        return 'None'


def get_crossmap_df(raw_df: pd.DataFrame, query_g: str, target_g: str, work_dir: str):
    print('get it')
    chain_fp = search_file(work_dir, 'chain.gz',
                           query_g=query_g, target_g=target_g, type='single')
    print(f"find chain file: {chain_fp}")
    global bx_ivl
    bx_ivl = chain.read_chain_file(chain_fp)[0]
    print('readed chain file')
    q_bed_fp = search_file(work_dir, 'bed', single_g=query_g, type='single')
    print(f"finded query_bed:{q_bed_fp}")
    t_bed_fp = search_file(work_dir, 'bed', single_g=target_g, type='single')
    print(f"finded target_bed:{t_bed_fp}")
    q_bed_df = pd.read_csv(q_bed_fp, sep='\t', header=None)
    q_bed_df.columns = ['chrom', 'start', 'end', 'genename', 'junk', 'strand']
    raw_bed_df = pd.merge(raw_df, q_bed_df, on='genename',
                          how='left').drop(columns=['junk'])
    global target_bed
    target_bed = pybedtools.BedTool(t_bed_fp)
    raw_bed_df['region'] = raw_bed_df.progress_apply(
        get_region, axis=1, args=(0.85,))
    raw_bed_df['crossmap'] = raw_bed_df['region'].progress_apply(
        string_bed_ivl, args=(0.8,))
    return raw_bed_df[['genename', 'crossmap']]
