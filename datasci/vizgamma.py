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


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('alpha', type=float)
    parser.add_argument('beta', type=float)
    parser.add_argument('outfile', nargs='?')
    parser.add_argument('-l', '--loc', default=0, type=float)

    args = parser.parse_args()

    if not args.outfile:
        # Use a graphical interface
        plt.switch_backend('TkAgg')

    alpha = args.alpha
    beta = args.beta
    loc = args.loc

    lines = plot_gamma_density(alpha, beta, loc=loc)
    ax = plt.gca()
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.text(xmax, ymax, 'E = %g' % expectation_gamma_beta(alpha, beta, loc),
            va='top', ha='right')
    ax.set_ylabel(r'Gamma(\alpha\,=%g; \beta\,=%g)' % (alpha, beta))
    plt.show(block=True)



if __name__ == '__main__':
    main()

