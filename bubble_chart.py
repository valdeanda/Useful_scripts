# coding: utf-8
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:        bubble_chart
# Purpose:     Program to creat a super easy bubble chart from a data frame
#
#
# @uthor:      acph - dragopoot@gmail.com
# @mods :      vda  - valdeanda@utexas.edu
# Created:     Jan, 2019
# Copyright:   (c) acph 2017
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------


#Libraries
import argparse
from argparse import RawDescriptionHelpFormatter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



epilog = """Example:
$  python3 bubble_chart.py data.tab """

parser = argparse.ArgumentParser(description=__doc__, epilog=epilog,
                                 formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('filename',
                    help="Input file dataframe i.e  abundances profile")
parser.add_argument('-im_format', '-f', default='pdf', type=str,
                    choices=['png', 'pdf', 'ps', 'eps', 'svg', 'tif', 'jpg'],
                    help='''Output format for images [pdf].''')
parser.add_argument('--im_res', '-r', default=300, type=int,
                    help='''Output resolution for images in
                    dot per inch (dpi) [dpi].''',
                    metavar='dpi')


args = parser.parse_args()
# END options #

df=pd.read_table(args.filename, index_col=0)
df=df.fillna(0)

df=df*100



def bubble_super_mega_and_simpe_plot(df, size_factor=25, xlabel='Columns',
                                     cmap='copper_r',alpha=0.3,
                                     ylabel='Rows'):
    #"Example of bubble plot, this function only works with small data, 0-10 approx,
    #if color map is None, scatter use simple color"
    n, m = df.shape
    X = df.values
    # here the processing for xc
    X_ = X * size_factor
    alpha = X_/X_.max()
    # now the scatters plots, one for rows
    xs = range(m)
    for y in range(n):
        ys = [y] * m
        if cmap:
            plt.scatter(xs, ys,
                        c=alpha[y],
                        # c='steelblue',
                        cmap=cmap,
                        vmin= alpha.min(),
                        vmax= alpha.max(),
                        s=X_[y])
        else:
            plt.scatter(xs, ys,
                        c='steelblue',
                        s=X_[y],
                       vmin=X_.min(),
                       vmax=X_.max())
    plt.xticks(np.arange(m), df.columns, rotation=90)
    plt.yticks(np.arange(n), df.index)
    plt.grid(linestyle='-', color='white')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.plot()
    if cmap:
        plt.colorbar()
    plt.tight_layout()
    plt.savefig(args.filename +"_bubbleplot." + args.im_format, dpi=args.im_res,bbox_inches='tight')


sns.set(font_scale=1)
sns.set_style("darkgrid")
plt.figure(figsize=(21,12))
plt.tight_layout()
bubble_super_mega_and_simpe_plot(df, 20, cmap='coolwarm_r', ylabel='Taxa', xlabel='Genes',alpha=0.05)

