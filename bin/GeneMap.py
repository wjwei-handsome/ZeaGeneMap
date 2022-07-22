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

import argparse

import sys
from typing import *
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from utils import utils
from utils.logger import pretty_print


LEGAL_EVIDENCES = ['rbh', 'crossmap', 'synteny', 'ortholog']

LOGO = """[bold green]
 _____                ___  ___            |
|  __ \               |  \/  |            | Author: Wenjie Wei
| |  \/ ___ _ __   ___| .  . | __ _ _ __  | Email: wjwei9908@gmail.com
| | __ / _ \ '_ \ / _ \ |\/| |/ _` | '_ \ | Contribute: Wenjie Wei;Songtao Gui;Han Xiao
| |_\ \  __/ | | |  __/ |  | | (_| | |_) || Department: HZAU-Maize-Bioinformatics
 \____/\___|_| |_|\___\_|  |_/\__,_| .__/ |
                                   | |    |
                                   |_|    |
demo: GeneMap.py -l list.txt -d work_dir -q GenomeA -t GenomeB -o output_prefix
"""


def get_args():
    # define arguments
    parser = argparse.ArgumentParser(
        description=None, prog='GeneMap.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required arguments
    required_parser = parser.add_argument_group('required arguments')
    required_parser.add_argument('-l', '--list', action='store',
                                 help='genename list file', required=True, dest='list')
    required_parser.add_argument('-d', '--dir', action='store',
                                 help='work dir', required=True, default='.', dest='dir')
    required_parser.add_argument('-q', '--query', action='store',
                                 help='query genome', required=True, dest='query_g')
    required_parser.add_argument('-t', '--target', action='store',
                                 help='target genome', required=True, dest='target_g')
    required_parser.add_argument('-o', '--output', action='store',
                                 help='output prefix', required=True, dest='output')

    # options arguments
    option_parser = parser.add_argument_group('optional arguments')
    option_parser.add_argument('-e', '--evidence',
                               help="evidence list", nargs='+', dest='evidence_l',
                               default=['rbh', 'crossmap', 'synteny', 'ortholog'])
    option_parser.add_argument('-r', '--ratio', action='store',
                               help='crossmap ratio', default=0.8, dest='ratio')
    option_parser.add_argument('-p', '--overlap', action='store',
                               help='overlap fraction', default=0.8, dest='overlap')

    return parser.parse_args()


def main(args):
    # judge the evidence
    utils.check_in_list(args.evidence_l, LEGAL_EVIDENCES,
                        print_list=True)

    # get raw genelist df
    raw_df = utils.get_raw_df(args.list)

    # add evidence and output raw df
    raw_evd_df = utils.get_raw_evd_df(args, raw_df, write=True)

    # process raw result
    utils.process_raw(args, raw_evd_df)


# TODO: README 4) 学会单元测试 5) 学会setup 6)发布

if __name__ == '__main__':
    pretty_print(LOGO)
    # get args
    args = get_args()
    main(args)
