#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pytest
from datasci.graphs import plottree, plottree_label_nodes, plottree_set_xlim
import ete3
import numpy as np
import pandas as pd
from dendro.any import ete3 as ete3_methods
import matplotlib as mpl
import matplotlib.pyplot as plt
import logging

get_items = ete3_methods.get_items
get_label = ete3_methods.get_label

logging.getLogger('datasci.graphs').setLevel(logging.DEBUG)


def get_answer(yesno_question):
    answer = None
    while answer not in ('n', 'y'):
        answer = input(yesno_question + ' [y/n] ').rstrip().lower()
    return answer


def test_triangles():
    #NOTE: the current implementation only shows a triangle if the collapsed node is a single child.
    tree = ete3.Tree('(((a)a0,(b)b0)x,(((c)c0,d)y,e)z)r;', format=1)
    collapsed = [(tree&name) for name in list('abc')]
    print('collapsed =', collapsed)
    lines, _, _ = plottree(tree, get_items, get_label, invert=False,
                           label_nodes=True, collapsed=collapsed)
    plt.show()
    plt.close()
    assert get_answer('Were triangles displayed?') == 'y'


class Test_edge_colormap:
    
    # NOTE
    # See thesis figures:
    # I.7 = ~/ws7/DUPLI_data93/alignments_analysis/fig/PhylTree-Simii_brlen-dS_edge-color-dSrate_fsahmmc_thesis.pdf
    #     computed in ~/ws7/DUPLI_data93/alignments_analysis/subtrees_stats/NBcells/Comparedating_make-control-ages-table.py
    # III.1 = ~/ws7/lateralite/data93/MMP21_presence-absence.pdf
    #     computed in ~/ws7/lateralite/find_candidates-Genomicus035-20200618.py

    tree = ete3.Tree('(((a)0,(b)1)x,(((c)2,d)y,e)z)r;', format=1)
    bin_values = pd.Series([1,1,0,0,1,1,0,1,0,1,0,1], index=list('a0b1c2dexyzr'), dtype=int)
    cat_values = pd.Series([1,1,0,0,-1,-1,1,0,-1,1,0,-1], index=list('a0b1c2dexyzr'), dtype=int)
    cont_values = pd.Series(np.linspace(-1, 1, 12), index=list('a0b1c2dexyzr'))

    @pytest.mark.parametrize('to_collapse', [None, list('abc')])
    def test_binary_edgevalues(self, to_collapse):
        collapsed = None
        if to_collapse is not None:
            collapsed = [(self.tree&name) for name in to_collapse]
        edge_cmap = plt.get_cmap('Paired', 2)
        lines, _, _ = plottree(self.tree, get_items, get_label,
                 invert=False,
                 edge_colors=self.bin_values,
                 edge_cmap=edge_cmap, label_nodes=True, collapsed=collapsed);
        lines.axes.legend([mpl.lines.Line2D([], [], color=edge_cmap(1)),
                           mpl.lines.Line2D([], [], color=edge_cmap(0))],
                         ['1', '0'], loc='upper left')
        plt.show()
        plt.close()
        assert get_answer('Were the colors correct?') == 'y'

    @pytest.mark.parametrize('to_collapse', [None, list('abc')])
    def test_categorical_edgevalues(self, to_collapse):
        collapsed = None
        if to_collapse is not None:
            collapsed = [(self.tree&name) for name in to_collapse]
        print('collapsed =', collapsed)
        fig, ax0 = plt.subplots()
        uniq_values = list(sorted(self.cat_values.unique()))
        n_values = len(uniq_values)
        edge_cmap = plt.get_cmap('Paired', n_values)
        edge_cmap.set_bad('k')
        lines, _, _ = plottree(self.tree, get_items, get_label,
                 invert=False,
                 edge_colors=self.cat_values,
                 edge_cmap=edge_cmap, label_nodes=True,
                 ax=ax0, collapsed=collapsed)
        edge_norm = mpl.colors.Normalize(self.cat_values.min(), self.cat_values.max())
        handles, labels = [], []
        for val in uniq_values:
            handles.append(mpl.lines.Line2D([], [], color=edge_cmap(edge_norm(val))))
            labels.append('%d' % val)
        lines.axes.legend(handles, labels, loc='upper left', title='Expected colors')

        #ax0.set_position(ax0.get_position().shrunk(0.5, 1.))
        #ax1 = fig.add_subplot(122)
        #ax1.imshow(np.array(uniq_values)[:,None], cmap=edge_cmap, norm=edge_norm)
        #ax1.set_yticks(np.arange(n_values))
        #ax1.set_yticklabels(uniq_values)
        #ax1.set_title('Expected colors')
        plt.show()
        plt.close()
        assert get_answer('Were the colors correct?') == 'y'


    @pytest.mark.parametrize('to_collapse', [None, list('abc')])
    def test_continuous_edgevalues(self, to_collapse):
        collapsed = None
        if to_collapse is not None:
            collapsed = [(self.tree&name) for name in to_collapse]
        edge_cmap = plt.get_cmap('viridis')
        edge_cmap.set_bad('grey')
        val_range = self.cont_values.max() - self.cont_values.min()
        edge_norm = mpl.colors.Normalize(self.cont_values.min() - 0.1*val_range,
                                         self.cont_values.max() + 0.1*val_range)

        lines,_,_ = plottree(self.tree, get_items, get_label,
                 invert=False,
                 label_nodes=True,
                 edge_colors=self.cont_values,
                 edge_cmap=edge_cmap, collapsed=collapsed)
        fig = plt.gcf()
        ax = plt.gca()
        lines.set_norm(edge_norm)
        cax,kw = mpl.colorbar.make_axes(ax, location='left', shrink=0.5, anchor=(1,1))
        cbar = fig.colorbar(lines, cax=cax, label='dS/My', **kw)
        plt.show()
        plt.close()
        assert get_answer('Were the colors correct?') == 'y'
