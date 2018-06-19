
import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import matplotlib.pyplot as plt
import numpy as np
import pylab
import warnings
warnings.filterwarnings("ignore")



wd = '/Users/songweizhi/Desktop'
ref_in = 'combined_mNC.fna'
os.chdir(wd)

for each in open('combined_ANI.txt'):
    if each.startswith('1'):

        print(each.strip())




