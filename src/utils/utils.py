from subprocess import Popen
import urllib
import gzip
import bz2
import sys
import os
import glob
import time
import pandas as pd
from tqdm import tqdm
from src.utils.logger import logger

tqdm.pandas(desc='Processing', colour='green')


def _warn(msg):
    logger.warning(f"{msg}")


def _error(msg):
    logger.error(f"{msg}")
    sys.exit(1)


def _info(msg):
    logger.info(f"{msg}")


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
            _error(
                f"{forward_fp} or {reverse_fp} not found or duplicated in your work dir : {abs_dir_path}")
    elif type == 'single':
        if query_g and target_g:
            chain_fp = f"{query_g}.{target_g}.{surfix}"
            chain_search = glob.glob(f"{abs_dir_path}/**/{chain_fp}")
            if len(chain_search) == 1:
                return chain_search[0]
            else:
                _error(
                    f"{chain_fp} not found or duplicated in your work dir : {abs_dir_path}")
        single_fp = f"{single_g}.{surfix}"
        single_search = glob.glob(f"{abs_dir_path}/**/{single_fp}")
        if len(single_search) == 1:
            return single_search[0]
        else:
            _error(
                f"{single_fp} not found or duplicated in your work dir : {abs_dir_path}")
    elif type == 'ortholog':
        ortholog_fp = f"rthogroups.{surfix}"
        ortholog_search = glob.glob(f"{abs_dir_path}/**/?{ortholog_fp}")
        if len(ortholog_search) == 1:
            return ortholog_search[0]
        else:
            _error(
                f"O[o]{ortholog_fp} not found or duplicated in your work dir : {abs_dir_path}")
    else:
        _error('type error, contact the author')


def nopen(f, mode="rb"):
    if not isinstance(f, str):
        return f
    if f.startswith("|"):
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
        if mode[0] == "r":
            return p.stdout
        return p
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
        else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
        else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
        else urllib.urlopen(f) if f.startswith(("http://", "https://", "ftp://")) \
        else open(f, mode)


def reader(fname):
    for l in nopen(fname):
        yield l.decode('utf8').strip().replace("\r", "")


def check_in_list(choice_list, legal_list, print_list: bool = True):
    if print_list:
        _info(f"choiced evidence: {choice_list}")
    if any(item not in legal_list for item in choice_list):
        _error(f'evidence error: please select in {legal_list}')


def get_time(func):
    """ a decorator to get the time of a function and return it """
    def wrapper(*args, **kwargs):
        start = time.time()
        res = func(*args, **kwargs)
        end = time.time()
        _info(f"Combine {func.__name__} evidenve cost {end - start}s ✅")
        return res
    return wrapper


def get_raw_df(genelist_fp: str) -> pd.DataFrame:
    # read input gene list
    try:
        raw_df = pd.read_csv(genelist_fp, sep='\t', header=None)
        _info(f"Read input genelist: {genelist_fp} ✅")
    except Exception as e:
        _error(f"Read {genelist_fp} failed, please check the file")
    if not raw_df.shape[1] == 1:
        _error('genelist file format error, except one column only')
    raw_df.columns = ['genename']
    return raw_df


@get_time
def rbh(args, raw_df):
    from src.blast import blast
    rbh_df = blast.get_rbh_df(args.query_g, args.target_g, args.dir)
    raw_rbh_merge = pd.merge(
        raw_df, rbh_df, on='genename', how='left')
    return raw_rbh_merge


@get_time
def ortholog(args, raw_df):
    from src.ortholog import ortholog
    ortholog_df = ortholog.get_ortholog_df(
        args.query_g, args.target_g, args.dir)
    raw_ortholog_merge = pd.merge(
        raw_df, ortholog_df, on='genename', how='left')
    if len(raw_ortholog_merge) != len(raw_df):
        _error("Ortholog evidence error")
    return raw_ortholog_merge


@get_time
def synteny(args, raw_df):
    from src.synteny import synteny
    rec_syn_df = synteny.get_rec_syn_df(
        args.query_g, args.target_g, args.dir)
    raw_syn_merge = pd.merge(
        raw_df, rec_syn_df, on='genename', how='left')
    return raw_syn_merge


@get_time
def crossmap(args, raw_df):
    from src.chain import crossmap
    crossmap_df = crossmap.get_crossmap_df(
        raw_df, args.query_g,
        args.target_g, args.dir,
        ratio=args.ratio, overlap=args.overlap
    )
    raw_crossmap_merge = pd.merge(
        raw_df, crossmap_df, on='genename', how='left')
    return raw_crossmap_merge


def get_raw_evd_df(args, raw_df, write: bool = True) -> pd.DataFrame:
    # get raw df
    # 1 rbh
    if 'rbh' in args.evidence_l:
        output_df = rbh(args, raw_df)
    else:
        output_df = raw_df
        _warn('rbh not in evidence list')
    # 2 ortholog
    if 'ortholog' in args.evidence_l:
        output_df = ortholog(args, output_df)
    else:
        _warn('ortholog not in evidence list')
    # 3 synteny
    if 'synteny' in args.evidence_l:
        output_df = synteny(args, output_df)
    else:
        _warn('synteny not in evidence list')
    # 4 crossmap
    if 'crossmap' in args.evidence_l:
        output_df = crossmap(args, output_df)
    else:
        _warn('crossmap not in evidence list')
    output_df.fillna('None', inplace=True)
    if write:
        _info(f'Writing raw result to {args.output}_raw.tsv')
        output_df.to_csv(
            f"{args.output}_raw.tsv", sep='\t', index=False, header=True)
    return output_df


def map_to_list(x):
    if '@' in x:
        lst = []
        for i in x.split(','):
            g = x.split('@')[0]
            if g not in lst:
                lst.append(g)
        return lst
    if x == 'NoneRegion':
        return ['None']
    else:
        return [item.strip() for item in x.split(',')]


def process_raw(args, raw_evd_df: pd.DataFrame) -> None:
    _info(f'Process raw result to final result')
    process_df = raw_evd_df.melt(
        id_vars=['genename'], var_name='evd', value_name='target')
    process_df['target'] = process_df['target'].progress_apply(map_to_list)
    process_df = process_df.explode('target')
    process_df = process_df.drop(
        process_df[process_df.target == 'None'].index).reset_index(drop=True)
    # process_df = process_df.drop(
    #     process_df[process_df.target == 'NoneRegion'].index).reset_index(drop=True)
    result_df = process_df.groupby(['genename', 'target'], as_index=False).agg({
        'evd': lambda x: ','.join(x)})
    _info(f'Writing final result to {args.output}_final.tsv')
    result_df.to_csv(f'{args.output}_final.tsv',
                     sep='\t', index=False, header=True)
