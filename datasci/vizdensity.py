#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Visualize the density of a gamma distribution with the provided parameters"""


import argparse as ap
import matplotlib as mpl
mpl.use('Agg')  # Do not connect to X11 by default
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st


def expectation_gamma(shape, scale, loc=0):
    return loc + shape*scale

def expectation_gamma_beta(shape, beta, loc=0):
    return loc + float(shape)/beta

def gamma_density(alpha=2, beta=0.5, scale=None, loc=0, ndots=250):
    if scale is None:
        # otherwise Beta is ignored
        scale = 1./beta

    # Go to 98%
    xmax = st.gamma.isf(0.02, alpha, scale=scale, loc=loc) #st.gamma.ppf(0.98, alpha)
    x = np.linspace(loc, xmax, ndots)

    return x, st.gamma.pdf(x, alpha, scale=scale, loc=loc)


def plot_gamma_density(alpha=2, beta=0.5, scale=None, loc=0, ndots=250):
    x, y = gamma_density(alpha, beta, scale, loc, ndots)
    return plt.plot(x, y, '-')


def print_gamma_density(alpha=2, beta=0.5, scale=None, loc=0):
    raise NotImplementedError
    #ndots =  # Measure number of columns
    x, y = gamma_density(alpha, beta, scale, loc, ndots)


def density(distrib, ndots=250):
    """Accepts a frozen distribution only"""
    xmin, xmax = distrib.a, distrib.b
    # Go to 98%
    if np.isinf(distrib.b):
        xmax = distrib.isf(0.02)
    if np.isinf(distrib.a):
        xmin = distrib.ppf(0.02)
    x = np.linspace(xmin, xmax, ndots)

    return x, distrib.pdf(x)


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parent = ap.ArgumentParser(add_help=False)
    parent.add_argument('-o', '--outfile')
    parent.add_argument('-l', '--loc', default=0, type=float)

    subp = parser.add_subparsers(title='command', dest='command')

    helpp = subp.add_parser('help', help='print the distrib docstring from scipy.stat')
    helpp.add_argument('distribname')

    gammap = subp.add_parser('gamma', parents=[parent])
    gammap.add_argument('alpha', type=float)
    gammap.add_argument('beta', type=float)

    generalp = subp.add_parser('general', parents=[parent])
    generalp.add_argument('distribname')
    generalp.add_argument('params', nargs='*')
    generalp.add_argument('-s', '--shape', default=1, type=float)

    args = parser.parse_args()

    if args.command == 'help':
        print(getattr(st, args.distribname).__doc__)
        return
    elif args.command == 'gamma':
        alpha = args.alpha
        beta = args.beta
        distribname = 'gamma'
        params = [args.alpha, args.loc, 1. / args.beta]
    else:
        distribname = args.distribname
        params = [float(x) for x in args.params]
        params += [args.loc, args.shape]

    # Freeze to specific parametrisation
    distrib = getattr(st, distribname)(*params)

    if not args.outfile:
        # Use a graphical interface
        plt.switch_backend('TkAgg')

    lines = plt.plot(*density(distrib), '-')
    ax = plt.gca()
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    quantiles = np.concatenate((distrib.ppf([0.05, 0.25, 0.5]),
                                distrib.isf([0.25, 0.05])))
    boxwidth = (ymax - ymin) / 10.
    ypad = (ymax-ymin) / 100.
    boxstats = dict(zip(('whislo', 'q1', 'med', 'q3', 'whishi'), quantiles))
    boxstats['mean'] = distrib.mean()
    box_elements = ax.bxp([boxstats], [-boxwidth / 2], widths=[boxwidth],
                          showfliers=False, showmeans=True,
                          vert=False, zorder=-1, manage_ticks=False)
    for quanti in quantiles:
        ax.vlines(quanti, 0+ypad, distrib.pdf(quanti)-ypad, colors=['.75'], linestyles='dashed')
    ticklabels =['q%d%%\n%.3g' % (percent, quanti) for percent, quanti in zip((5, 25, 50, 75, 95), quantiles)]
    ax.set_xticks(quantiles)
    ax.set_xticklabels(ticklabels, fontsize='small')

    ax.legend(box_elements['means'], ['mean = %g' % distrib.mean()])

    if args.command == 'gamma':
        ax.set_ylabel(r'$\Gamma(\alpha\,=%g; \beta\,=%g)$' % (alpha, beta))
    else:
        ax.set_ylabel(r'$%s(%s, s=%g, loc=%g)$' % (distribname, ', '.join(args.params), args.shape, args.loc))

    # Styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.show(block=True)


if __name__ == '__main__':
    main()
