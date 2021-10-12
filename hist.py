# fie coding: utf-8
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     hist.py
# Purpose:  plot depth file histogram of big sequencing data
# Authors:  vydat - valdeanda@ciencias.unam.mx
# Created:  2018
# Licence:  GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
import argparse
from argparse import RawDescriptionHelpFormatter
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns

epilog = """Example:
$  python3 histplot.py  depthfile.tsv"""

parser = argparse.ArgumentParser(description=__doc__, epilog=epilog,
                                 formatter_class=RawDescriptionHelpFormatter)

parser.add_argument('filename',
                    help="lenght file")

parser.add_argument('-im_format', '-f', default='pdf', type=str,
                    choices=['png', 'pdf', 'ps', 'eps', 'svg', 'tif',
                             'jpg'],
                    help='''Output format for images [pdf].''')

parser.add_argument('--im_res', '-r', default=300, type=int,
                    help='''Output resolution for images in
                    dot per inch (dpi) [dpi].''',
                    metavar='dpi')

args = parser.parse_args()

filename = args.filename

df =pd.read_csv(filename,sep='\t',index_col=0)

sns.set(font_scale=1.5)
fig = plt.figure(figsize=(12,7))
ax = sns.distplot(df, rug=True, rug_kws={"color": "b"},
                  kde_kws={"color": "k", "lw": 3, "label": "KDE"},
                  hist_kws={"histtype": "step", "linewidth": 3,
                  "alpha": 1, "color": "b"})
plt.title(str(args.filename), fontweight='bold')
plt.ylabel('Density', fontweight='bold')
plt.xlabel('Length', fontweight='bold')
plt.tight_layout()
plt.savefig(args.filename + "_hisplot." + args.im_format, dpi=args.im_res, bbox_inches='tight')
