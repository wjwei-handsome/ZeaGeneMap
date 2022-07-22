from dataclasses import dataclass, field
from typing import *
import sys
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from bin.GeneMap import main


@dataclass
class TestArgs:
    list: str = 'tests/data/Zm-B73v5.genelist'
    dir: str = '.'
    query_g: str = 'Zm-B73v5'
    target_g: str = 'Zv-TEO01'
    output: str = 'tests/data/test'
    evidence_l: List = field(default_factory=lambda: [
                             'rbh', 'crossmap', 'synteny', 'ortholog'])
    ratio: float = 0.8
    overlap: float = 0.8


TestArgs = TestArgs()

main(TestArgs)
