from subprocess import Popen
import urllib
import gzip
import bz2
import sys
import logging
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
        if query_g and target_g:
            chain_fp = f"{query_g}.{target_g}.{surfix}"
            print(chain_fp)
            chain_search = glob.glob(f"{abs_dir_path}/**/{chain_fp}")
            print(chain_search)
            if len(chain_search) == 1:
                return chain_search[0]
            else:
                print('ss')
                exit(1)
        single_fp = f"{single_g}.{surfix}"
        single_search = glob.glob(f"{abs_dir_path}/**/{single_fp}")
        if len(single_search) == 1:
            return single_search[0]
        else:
            exit(1)
    elif type == 'ortholog':
        search = glob.glob(f"{abs_dir_path}/**/?rthogroups.{surfix}")
        if len(search) == 1:
            return search[0]
        else:
            exit(1)
    else:
        logging.error('type error')


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
