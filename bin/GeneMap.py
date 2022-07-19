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

import logging
import argparse


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
print(args.evidence_l)
