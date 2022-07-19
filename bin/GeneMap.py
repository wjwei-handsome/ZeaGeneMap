#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :GeneMap.py
@说明        :
@时间        :2022/07/18 21:02:56
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''

# from blast.blast import get_rbh_result
import argparse
import logging
import sys
import pandas as pd
from types import *
import os
from tqdm import tqdm
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)


tqdm.pandas()


def get_args():
    parser = argparse.ArgumentParser(
        description=__doc__, prog='GeneMap.py')
    parser.add_argument('-l', '--list', action='store',
                        help='genename list file', required=True, dest='list')
    parser.add_argument('-d', '--dir', action='store',
                        help='work dir', required=True, default='.', dest='dir')
    parser.add_argument('-q', '--query', action='store',
                        help='query genome', required=True, dest='query_g')
    parser.add_argument('-t', '--target', action='store',
                        help='target genome', required=True, dest='target_g')
    parser.add_argument('-e', '--evidence',
                        help='evidence list', nargs='+', dest='evidence_l',
                        default=['rbh', 'crossmap', 'synteny', 'ortholog'])
    # parser.add_argument('-d', '--dir', action='store_true',
    #                     help='work dir', required=True)
    return parser.parse_args()


args = get_args()

# judge the evidence
evidence_list = args.evidence_l
print(evidence_list)
if any(e not in ['rbh', 'crossmap', 'synteny', 'ortholog'] for e in evidence_list):
    logging.error('evidence error')
    exit(1)

# write result
genelist_fp = args.list

raw_df = pd.read_csv(genelist_fp, sep='\t', header=None)
assert raw_df.shape[1] == 1, 'genelist file format error'
raw_df.columns = ['genename']


if 'rbh' in args.evidence_l:
    from src.blast import blast
    rbh_df = blast.get_rbh_df(args.query_g, args.target_g, args.dir)
    raw_df['rbh'] = pd.merge(
        raw_df, rbh_df, on='genename', how='left')['target']
    raw_df.fillna('None', inplace=True)
    raw_df.to_csv('test_out.tsv', sep='\t', index=False, header=False)
else:
    print('rbh not in evidence list')
if 'ortholog' in args.evidence_l:
    ...
else:
    print('ortholog not in evidence list')
if 'synteny' in args.evidence_l:
    ...
else:
    print('synteny not in evidence list')
if 'crossmap' in args.evidence_l:
    ...
else:
    print('crossmap not in evidence list')
