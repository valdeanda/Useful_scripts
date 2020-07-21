#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:
# Purpose:
#
# @uthor:   acph - dragopoot@gmail.com
#
# Created:
# Copyright:   (c) acph 2015
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
""" This script create a cluster map """

import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.cluster.hierarchy as sch



def argparser():
    """Arguments"""
    epilog = """Example:
    $  python3 heatmap.py data.heatmap.tsv """

    parser = argparse.ArgumentParser(description=__doc__, epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename',
                        help="Input file derived mebs_output with classification")
    parser.add_argument('-f', '--im_format', default='png', type=str,
                        choices=['png', 'pdf', 'ps', 'eps', 'svg', 'tif', 'jpg'],
                        help='''Output format for images [png].''')
    parser.add_argument('-r', '--im_res', default=300, type=int,
                        help='''Output resolution for images in
                        dot per inch (dpi) [dpi].''',
                        metavar='dpi')
    # parser.add_argument('-o', '--outfig') # outfile
    # parser.add_argument('') # cluster files
    # parser.add_argument('') # distance metric
    # parser.add_argument('') # linkage method
    # parser.add_argument('') # columns labels
    # parser.add_argument('') # colormap
    args = parser.parse_args()
    return args


def main():
    """Main function"""
    args = argparser()
    matplotlib.rcParams['lines.linewidth'] = 0.4
    # load data
    print('[INFO] Loading data')
    data = pd.read_csv(args.filename,
                       sep='\t',
                       index_col=[0, 1])

    # *************************************************************************
    # *                         data scaling normaliation                     *
    # *************************************************************************
    # to do
    # data_scaling = (data.T - data.T.min())/(data.T.max() - data.T.min())
    # data = data_scaling.T
    # *************************************************************************
    # *               Fix columns (rows) with 0 or show error and             *
    # *                            warning to the user                        *
    # *************************************************************************
    # Fix for 0 columns
    print("[WARN] Fixing all 0s columns")
    data = data.loc[:, data.sum() != 0]

    # *************************************************************************
    # *                   Calculate distnaces - add to options                *
    # *************************************************************************
    # Rows distances
    d_metric = 'braycurtis'
    linkage_m = 'average'
    metadist = sch.distance.pdist(data, metric=d_metric)
    metalink = sch.linkage(metadist, method=linkage_m)
    metalink = metalink.clip(0, metalink.max() + 1)
    # columns distances
    profdist = sch.distance.pdist(data.T, metric=d_metric)
    proflink = sch.linkage(profdist, method=linkage_m)
    proflink = proflink.clip(0, proflink.max() + 1)

    ############
    # Plotting #
    ############
    print('[INFO] Plotting ...')
    # - Figure setup
    xf = 6.7
    yf = 8.6
    fig = plt.figure(figsize=(xf, yf))
    # Axes positions
    # # Axes without column names
    # posm = [0.01, 0.01, 0.2, 0.82]
    # posp = [0.24, 0.84, 0.67, 0.15]
    # posmat = [0.24, 0.01, 0.67, 0.82]
    # posm_colors = [0.215, 0.01, 0.02, 0.82]
    # poscbar = [0.92, 0.01, 0.015, 0.40]

    # new with labesl
    posm = [0.01, 0.23, 0.2, 0.62]
    posp = [0.24, 0.855, 0.67, 0.14]
    posmat = [0.24, 0.23, 0.67, 0.62]
    posm_colors = [0.215, 0.23, 0.02, 0.62]
    # poscbar = [0.94, 0.01, 0.02, 0.41]
    poscbar = [0.92, 0.23, 0.015, 0.30]
    poslegend = [0.01, 0.84, 0.23, 0.15]

    # colors for dendograms
    sch.set_link_color_palette(['#1f77b4', '#ff7f0e',
                                '#2ca02c', '#d62728',
                                '#9467bd', '#8c564b',
                                '#e377c2', '#7f7f7f',
                                '#bcbd22', '#17becf'])

    # # - rows dendogram
    meta_ax = fig.add_axes(posm, frameon=False)
    metadend = sch.dendrogram(metalink,
                              color_threshold=0.2 * max(metalink[:, 2]),
                              orientation='left')
    meta_ax.set_xticks([])
    meta_ax.set_yticks([])

    # # - columns dendogram
    prof_ax = fig.add_axes(posp, frameon=False)
    profdend = sch.dendrogram(proflink,
                              color_threshold=0.2 * max(proflink[:, 2]),
                              orientation='top')
    prof_ax.set_xticks([])
    prof_ax.set_yticks([])

    # # - Matrix - HEATMAP
    matrix_ax = fig.add_axes(posmat)
    mat = data.get_values()
    mmask = metadend['leaves']
    pmask = profdend['leaves']
    mat = mat[mmask, :]
    mat = mat[:, pmask]
    im = matrix_ax.matshow(mat, aspect='auto', origin='lower', cmap='viridis')
# ****************************************************************************
# *  Here we can add options to show labels - needs to modify axes positions *
# ****************************************************************************
    # Etiquetas
    matrix_ax.set_yticks([])
    matrix_ax.set_xticks(range(len(data.columns)))
    matrix_ax.set_xticklabels(data.columns[pmask], rotation=90,
                              fontsize=5, color='k')
    matrix_ax.xaxis.set_ticks_position('bottom')

    # # - Colorbar
    colorbar_ax = fig.add_axes(poscbar)
    cb = plt.colorbar(im, cax=colorbar_ax)
    cb.set_label('Completeness', fontsize='x-small')
    cb.ax.tick_params(labelsize='xx-small')

    # # - Color code
    # Get colors from the first element in the index
    general_index = sorted(pd.MultiIndex.to_frame(data.index)[0].unique(),
        key=lambda x: int(x[1:]))
    color_as = {}
    for i, v in enumerate(range(len(general_index))):
        icolor = plt.cm.tab20(v / len(general_index))
        color_as[general_index[i]] = icolor
    # color vector
    color_vec = [color_as[i[0]] for i in data.index]
    color_vec = np.array(color_vec)
    color_vec = color_vec[mmask]
    color_ax = fig.add_axes(posm_colors, frameon=False)
    lefts = range(0, len(color_vec), 1)
    height = np.ones(len(color_vec))
    width = 1
    metabars = color_ax.barh(lefts, height, width, color=color_vec,
                             edgecolor=color_vec)
    # Can you use matshow, pcolor or imshow?
    # im_col = color_ax.matshow(color_mat, aspect='auto',
    #                           origin="lower")
    # color_ax.set_xlim(-0.5, 0.5)
    color_ax.set_xticks([])
    color_ax.set_yticks([])
    color_ax.set_ylim((0, len(color_vec)))

    # # - Legend
    legend_ax = fig.add_axes(poslegend, frameon=False)
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])
    patches = []
    for name, color_ in color_as.items():
        p = mpatches.Patch(color=color_, label=name)
        patches.append(p)
    plt.legend(handles=patches, fancybox=True, fontsize='xx-small',
               loc=2, framealpha=0.75)
    # # - show
    # plt.show()
    # # - Save Figure
    figname = 'heatmap.{}'.format(args.im_format)
    fig.savefig(figname, dpi=args.im_res)

# ****************************************************************************
# *                               Save cluster data                          *
# ****************************************************************************
    # # - Save cluster data!!!!
    # clustdata = data.iloc[mmask, pmask]
    # clustdata.to_csv('matdata_clust_pru.txt', sep='\t')



if __name__ == '__main__':
    main()
    print('[INFO] End!!!. Thanks for using :D')




# ****************************************************************************
# *                                T R A S H                                 *
# ****************************************************************************




    # # loadign and creating top pfam
    # def labels_names(data, namelist):
    #     """Plot in x label the names in list

    #     Arguments:
    #     - `data`: panda DataFrame
    #     - `namelist`: a python list of names. Must be contained in
    #     pandas data.columns
    #     """
    #     labels = []
    #     for col in data.columns:
    #         if col in namelist:
    #             labels.append(col)
    #         else:
    #             labels.append('')
    #     return np.array(labels)

    # top_pfam = pd.read_table('top_pfam_profiles.txt', index_col=0)
    # top_p = list(top_pfam.index)
    # labs = labels_names(data, top_p)
    # labs = labs[pmask]
