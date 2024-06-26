#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Various statistical and graph routines.

For DataFrame operations, see `.dataframe_recipees`.

Warning: some of them can only be run in a notebook (because using display_html
or Pandas dataframe styling.)"""

from io import StringIO
from functools import partial
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from IPython.display import display, display_html
import warnings

import scipy.stats as stats
import scipy.cluster.hierarchy as hclust

import statsmodels.api as sm
from sklearn.decomposition import PCA, FactorAnalysis
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score

from datasci.graphs import scatter_density, \
                           plot_cov, \
                           plot_features_radar, \
                           heatmap_cov
from datasci.stats import normal_fit, cov2cor, r_squared, adj_r_squared
from datasci.dataframe_recipees import centered_background_gradient, magnify, \
                                        matplotlib_stylebar
import logging

#logfmt = "%(levelname)-7s:%(name)s:%(funcName)-20s:%(message)s"
#logf = logging.Formatter(logfmt)
#
#try:
#    from UItools import colorlog
#    #clogfmt = "$LVL%(levelname)-7s:$RESET${white}l.%(lineno)3s:%(funcName)-20s:$RESET%(message)s"
#    clogfmt = "$LVL%(levelname)-7s:$RESET${white}l.%(lineno)2s:$RESET%(message)s"
#    colorlogf = colorlog.ColoredFormatter(clogfmt)
#except ImportError:
#    colorlogf = logf

logger = logging.getLogger(__name__)


# Variable transformation

def notransform(x):
    return x

#notransform.__name__ = "%s"
sqrt = np.sqrt
log = np.log10

def logneg(x):
    return -np.log10(-x)
logneg.__name__ = "-log10(-%s)"

def sqrtneg(x):
    return -sqrt(-x)

def logsqrt(x):
    return np.log10(sqrt(x))
logsqrt.__name__ = "log10(sqrt(%s))"

# These should be classes, with attributes mainfunc, inc, offset, neg, and a __call__ method.
def make_sqrtpostransform(inc=0):
    sqrtinc = lambda x: sqrt(x - x.min() + inc)
    sqrtinc.__name__ = "sqrt(%g-min+%%s)" % inc
    return sqrtinc

def make_logtransform_inc(inc=0):
    loginc = lambda x: np.log10(x + inc)
    loginc.__name__ = "log10(%g+%%s)" % inc
    return loginc

def make_lognegtransform_inc(inc=0):
    logneginc = lambda x: np.log10(-x + inc)
    #FIXME: should be
    #logneginc = lambda x: -np.log10(-x + inc)
    logneginc.__name__ = "log10(%g-%%s)" % inc
    return logneginc

def make_logpostransform_inc(inc=0):
    loginc = lambda x: np.log10(x - x.min() + inc)
    loginc.__name__ = "log10(%g-min+%%s)" % (inc)
    return loginc


def make_best_logtransform(var, quantile=0.05):
    inc = 0
    if (var < 0).any():
        if (var > 0).any():
            # Both negative and positive values
            current_logtransform = make_logpostransform_inc
        else:
            current_logtransform = make_lognegtransform_inc
            if (var==0).any():
                inc = -var[var!=0].quantile(1-quantile)
    else:
        current_logtransform = make_logtransform_inc
        if (var==0).any():
            inc = var[var!=0].quantile(quantile)

    return current_logtransform(inc)


def discretizer(bins=[0], right=True):
    """"""
    def discretize(x):
        return np.digitize(x, bins, right=right)
    if len(bins) == 1:
        discretize.__name__ = 'binarize_%g' % bins[0]
    return discretize


def binarize(x):
    return np.digitize(x, [0], right=True)

binarize0 = binarize

def binarize1(x):
    return np.digitize(x, [1], right=False)


def zscore(x, ddof=1):
    return (x - x.mean()) / x.std(ddof=ddof) # if x.std() else 1) )

def standardize(x, ddof=1):
    """Rescale a feature to have std=1. Use this instead of zscore for binary features."""
    return x/x.std(ddof=ddof)

def unregress(y, x):
    """Return the residuals of Y from the linear regression (with Intercept)
    against X"""
    fit = sm.OLS(y, sm.add_constant(x)).fit()
    logger.debug('unregress Nobs=%d: a=%g; b=%g', fit.nobs, fit.params['const'], fit.params[1])
    return fit.resid

def unregress0(y, x):
    """Return the residuals of Y from the linear regression (with Intercept 0)
    against X"""
    return sm.OLS(y, x).fit().resid

def decorrelator(decorrfunc, data, src_data=None, *args, **kwargs):
    """
    *args: tuples of ('y', 'x') of variables in src_data;
    **kwargs: dict(y='x') of variables in src_data;

    Apply `decorrfunc(src_data[y], src_data[x])` to the new data.
    """
    if src_data is None:
        src_data = data

    kwargs.update({'R%s' % v[0]: v for v in args})

    return data.assign(**{newvar: decorrfunc(src_data[var],src_data[corrvar])
                          for newvar, (var, corrvar) in kwargs.items()})\
               .drop([var for newvar, (var,_) in kwargs.items()
                      if newvar != var],
                     axis=1, errors='ignore')


def check_decorrelator(decorrfunc, data, src_data=None, *args, **kwargs):
    decorrelated = decorrelator(decorrfunc, data, src_data, *args, **kwargs)
    if src_data is None:
        src_data = data
    kwargs.update({'R%s' % v[0]: v for v in args})

    # Display the Old VS New correlation coefficient.
    for newvar, (var, corrvar) in kwargs.items():
        old_r, old_pval = stats.spearmanr(src_data[var], src_data[corrvar])
        new_r, new_pval = stats.spearmanr(decorrelated[newvar], src_data[corrvar])
        if abs(old_r) <= abs(new_r):
            logger.warning('Did not decrease correlation: %s / %s: '
                           'abs(Spearman R) %.4f -> %.4f', var, corrvar,
                           old_r, new_r)
    return decorrelated

# These should be classes
class Decorrelator(object):
    def __init__(self, decorrfunc):
        self.decorrfunc = decorrfunc
    def __call__(self, data, src_data=None, *args, **kwargs):
        if src_data is None:
            src_data = data

        kwargs.update({'R%s' % v[0]: v for v in args})

        return data.assign(**{newvar: self.decorrfunc(src_data[var],src_data[corrvar])
                              for newvar, (var, corrvar) in kwargs.items()})\
                   .drop([var for newvar, (var,_) in kwargs.items()
                          if newvar != var],
                         axis=1, errors='ignore')
    def __str__(self):
        return '%s(%s)' % (self.__class__.__name__,
                getattr(self.decorrfunc, '__name__', str(self.decorrfunc)))
    def __repr__(self):
        return '<%s at 0x%x>' % (str(self), id(self))

class checkDecorrelator(Decorrelator):
    def __call__(self, data, src_data=None, *args, **kwargs):
        decorrelated = super().__call__(data, src_data, *args, **kwargs)
        if src_data is None:
            src_data = data
        kwargs.update({'R%s' % v[0]: v for v in args})

        # Display the Old VS New correlation coefficient.
        for newvar, (var, corrvar) in kwargs.items():
            old_r, old_pval = stats.spearmanr(src_data[var], src_data[corrvar])
            new_r, new_pval = stats.spearmanr(decorrelated[newvar], src_data[corrvar])
            if abs(old_r) <= abs(new_r):
                logger.warning('Did not decrease correlation: %s / %s: '
                               'abs(Spearman R) %.4f -> %.4f', var, corrvar,
                               old_r, new_r)
        return decorrelated
    

decorrelate = partial(check_decorrelator, np.divide)        # = checkDecorrelator(np.divide)
decorrelatelogs = partial(check_decorrelator, np.subtract)  # = checkDecorrelator(np.subtract)

def zscore_divide(v, cv): return zscore(np.divide(v,cv))
def zscore_subtract(v, cv): return zscore(np.subtract(v,cv))
def standardize_subtract(v, cv): return standardize(np.subtract(v,cv))
def zscore_unregress(v, cv): return zscore(unregress(v,cv))
def zscore_unregress0(v, cv): return zscore(unregress0(v,cv))

renorm_decorrelate = partial(check_decorrelator, zscore_divide)  # = checkDecorrelator(zscore_divide)
renorm_decorrelatelogs = partial(check_decorrelator, zscore_subtract)
check_unregress = partial(check_decorrelator, unregress)
renorm_unregress = partial(check_decorrelator, zscore_unregress)
renorm_unregress0 = partial(check_decorrelator, zscore_unregress0)


def renorm_bestlog_and_divide(v1, v2):
    tr_v1 = make_best_logtransform(v1)
    tr_v2 = make_best_logtransform(v2)
    if not (tr_v1.__name__.startswith('log10(0+') and
            tr_v2.__name__.startswith('log10(0+')):
        logger.info('Transforms: %s, %s', tr_v2.__name__, tr_v2.__name__)
    return zscore(tr_v1(v1) / tr_v2(v2))


def renorm_bestlog_and_unregress(v1, v2):
    tr_v1 = make_best_logtransform(v1)
    tr_v2 = make_best_logtransform(v2)
    if not (tr_v1.__name__.startswith('log10(0+') and
            tr_v2.__name__.startswith('log10(0+')):
        logger.info('Transforms: %s, %s', tr_v2.__name__, tr_v2.__name__)
    return zscore(unregress(tr_v1(v1), tr_v2(v2)))


def make_retro_trans(transform):
    if transform == 'notransform':
        retro_trans = notransform
    elif transform == 'log10':
        def retro_trans(x): return 10**x
    elif transform == '-log10(-%s)':
        def retro_trans(x): return -10**(-x)
    elif transform.startswith('log10('):
        if transform.endswith('-min+%s)'):
            inc = float(transform.replace('log10(', '').replace('-min+%s)', ''))
            miny = self.alls[y].min()
            def retro_trans(x):
                return 10**x + miny + inc
        elif transform.endswith('+%s)'):
            inc = float(transform.replace('log10(', '').replace('+%s)', ''))
            def retro_trans(x): return 10**x - inc
        elif transform.endswith('-%s)'):
            inc = float(transform.replace('log10(', '').replace('-%s)', ''))
            def retro_trans(x): return -10**(x) - inc
    elif transform == 'sqrt':
        def retro_trans(x): return x*x
    elif transform == 'sqrtneg':
        def retro_trans(x): return -x*x
    elif transform.startswith('sqrt('):
        if transform.endswith('-min+%s)'):
            miny = self.alls[y].min()
            inc = float(transform.replace('sqrt(', '').replace('-min+%s)', ''))
            def retro_trans(x):
                return x*x + miny + inc
    else:
        logger.error('No reverse transformation known for %r' % transform)
        retro_trans = notransform
    return retro_trans


def test_transforms(alls, variables, figsize=(14, 5), out=None, widget=None):
    """:param: out: file-like object with `.write()` method (given to `print`).
    If out.output, out.html or out.show exist, they replace
    display, display_html and plt.show.
    """
    #fig, axes = plt.subplots(len(variables),3, figsize=(22, 5*len(variables)))
    #if out is not None:
        # Replace the functions: `print`, `display`, `plt.show`
        #if widget:
        #    raise ValueError("arguments `out` and `widget` are incompatible")
        #FIXME: Implementation could probably be nicer with something like Mock.patch.
        # From the datasci.savior.*Report classes
    disp = getattr(out, 'output', display)
    disp_html = getattr(out, 'html', display_html)
    show = getattr(out, 'show', plt.show)  # Bad things may happen with unexpected classes.

    nbins = 50

    suggested_transform = {}

    iter_features = enumerate(variables)
    if widget is not None:
        #iter_features = zip(iter_features, widget())  # Example: widget=lambda: generate_slideshow(hr)
        iter_features = widget(iter_features)  # Example: widget=slideshow_generator(hr)
        #try:
        #    from UItools.jupytertricks import make_tabs
        #    iter_features = make_tabs(iter_features, lambda item: item[1])
        #except ImportError as err:
        #    logger.error()
    for i, ft in iter_features:

        var = alls[ft]

        transform_skews = {}
        # Plot original variable distribution 
        if var.dtype != float:
            print("Variable %r not continuous: %s" % (ft, var.dtype), file=out)
            if not np.issubdtype(var.dtype, np.integer):
                if np.issubdtype(var.dtype,  bool):
                    if var.min() == var.max():
                        print("Boolean variable %r is constant. Skipping." % ft, file=out)
                    else:
                        print("Boolean variable %r interpreted as numeric" % ft, file=out)
                        suggested_transform[ft] = binarize
                else:
                    print("Variable %r does not seem to be numeric/bool. Skipping" % ft, file=out)
                continue
        if var.min() == var.max():
            print("Variable %r is constant. Skipping." % ft, file=out)
            # Why does it appear in `suggested_transform`?
            continue
        
        infinite_vals = np.isinf(var)  # Ignore NAs
        if infinite_vals.any():
            print("%d infinite values in %r (observations ignored)." % (infinite_vals.sum(), ft), file=out)
            var = var[~infinite_vals]

        fig, axes = plt.subplots(1, 3, figsize=figsize)
        #axes[i, 0].set_ylabel(ft)
        axes[0].set_ylabel(ft, fontsize='large')

        axes[0].hist(var, bins=nbins, density=True)
        axes[0].plot(*normal_fit(var), '-')
        axes[0].set_title("original")
        _, xmax0 = axes[0].get_xlim()
        _, ymax0 = axes[0].get_ylim()
        
        varskew = var.skew()
        text = "Skew: %g\nKurtosis: %g\n" % (varskew, var.kurt())
        transform_skews[notransform] = varskew

        current_logtransform = None
        current_signtransform = None  #TODO
        current_sqrttransform = None
        if (var < 0).any():
            if (var > 0).any():
                print("Variable %r has negative and positive values. Shifting to positive." % ft, file=out)
                text += "Negative and positive values. Shifting to positive.\n"
                var -= var.min()
                current_logtransform = make_logpostransform_inc()
                current_sqrttransform = make_sqrtpostransform()
            else:
                print("Variable %r converted to positive values" % ft, file=out)
                text += "Converted to positive values.\n"
                var = -var
                current_logtransform = logneg  # What if null values?
                current_sqrttransform = sqrtneg
        else:
            current_logtransform = log
            current_sqrttransform = sqrt
        
        # Plot log-transformed distribution
        with warnings.catch_warnings(record=True) as w:
            logtransformed_var = np.log10(var)
        if w:
            if (issubclass(w[-1].category, RuntimeWarning)
                and "divide by zero encountered in log10" in str(w[-1].message)):
                zero_vals = True
                n_zero_vals = (var == 0).sum()
            else:
                warnings.warn(w[-1].category(w[-1].message))

        # Not defined for log.
        n_infinite_vals = (~np.isfinite(logtransformed_var)).sum()
        logtransformed_var = logtransformed_var[np.isfinite(logtransformed_var)]
        
        logskew, logkurt = logtransformed_var.skew(), logtransformed_var.kurt()
        logtext = "Skew: %g\nKurtosis: %g\n" % (logskew, logkurt)

        axes[1].hist(logtransformed_var, bins=nbins, density=True, alpha=0.5)
        axes[1].plot(*normal_fit(logtransformed_var), '-', alpha=0.5)
        
        if not n_infinite_vals:
            transform_skews[current_logtransform] = logskew
        else:
            make_logtransform = make_logtransform_inc if current_logtransform is log \
                                else make_logpostransform_inc
            
            suggested_increment = logtransformed_var.quantile(0.05)
            print("%s: Nb of not finite values: %d. Suggested increment: 10^(%g) = %g"\
                  % (ft, n_infinite_vals, suggested_increment, 10**suggested_increment),
                  file=out)
            text += "%d not finite values. Suggested increment: 10^(%g)"\
                     % (n_infinite_vals, suggested_increment)

            logtransformed_inc_var = np.log10(var + 10**suggested_increment) #1.) 
            twin_ax1 = axes[1].twinx()
            twin_ax1.hist(logtransformed_inc_var, bins=nbins, density=True,
                          alpha=0.5, label="Incremented", color="#78a86a") # Green
            twin_ax1.plot(*normal_fit(logtransformed_inc_var), '-', alpha=0.5)
            
            logtransformed_inc1_var = np.log10(var + 1)
            twin2_ax1 = axes[1].twinx()
            twin2_ax1.hist(logtransformed_inc1_var, bins=nbins, density=True,
                           alpha=0.5, label="Increment+1", color="#a86a78") # Red
            twin2_ax1.plot(*normal_fit(logtransformed_inc1_var), '-', alpha=0.5)
            
            log_inc_varskew = logtransformed_inc_var.skew()
            logtext_inc = ("Skew     (+inc): %g\nKurtosis (+inc): %g\n"
                           % (log_inc_varskew, logtransformed_inc_var.kurt()))
            transform_skews[make_logtransform(10**suggested_increment)] = log_inc_varskew

            log_inc1_varskew = logtransformed_inc1_var.skew()
            logtext_inc1 = ("Skew     (+1): %g\nKurtosis (+1): %g\n"
                            % (log_inc1_varskew, logtransformed_inc1_var.kurt()))
            transform_skews[make_logtransform(1)] = log_inc1_varskew

        xmin1, xmax1 = axes[1].get_xlim()
        _, ymax1 = axes[1].get_ylim()

        axes[1].set_title("log10 transformed")
        
        sqrttransformed_var = np.sqrt(var)
        sqrttransformed_varskew = sqrttransformed_var.skew()
        sqrttext = "Skew: %g\nKurtosis: %g\n" % (sqrttransformed_varskew,
                                                 sqrttransformed_var.kurt())
        transform_skews[current_sqrttransform] = sqrttransformed_varskew

        axes[2].hist(sqrttransformed_var, bins=nbins, density=True)
        axes[2].plot(*normal_fit(sqrttransformed_var), '-')
        axes[2].set_title("Square root transformed")
        _, xmax2 = axes[2].get_xlim()
        _, ymax2 = axes[2].get_ylim()
        
        axes[0].text(xmax0, ymax0, text, va="top", ha="right")

        xpos1 = xmax1 if logskew>0 else xmin1
        ha1 = 'right' if logskew>0 else 'left'
        if n_infinite_vals:
            axes[1].text(xpos1, ymax1, logtext_inc, va="top", ha=ha1, color='#78a86a')
            axes[1].text(xpos1, 0.85*ymax1, logtext_inc1, va="top", ha=ha1, color='#a86a78')
        else:
            axes[1].text(xpos1, ymax1, logtext, va="top", ha=ha1)

        axes[2].text(xmax2, ymax2, sqrttext, va="top", ha="right")

        #fig.show(warn=False)
        show()
        plt.close(fig)
        suggested_transform[ft] = sorted(transform_skews.items(),
                                         key=lambda x: abs(x[1]))[0][0]

    return suggested_transform  #, figs


def randomforest_regression(X, y, cv=10, n_estimators=100, max_depth=5, out=None, **kwargs):
    RF = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth)
    RFfit = RF.fit(X, y)
    ft_importance = pd.Series(RFfit.feature_importances_,
                              index=X.columns).sort_values(ascending=False)
    print("R² =", RFfit.score(X, y), file=out)
    print(ft_importance, file=out)  # Does not give the direction of the effect.

    # a list of R² for each fit without a subsample
    RFcrossval = cross_val_score(RF, X, y, cv=cv, n_jobs=-1)
    print(RFcrossval, file=out)

    return RFfit


def print_pca_loadings(pca, features):
    """Works with sklearn.decomposition.PCA or FactorAnalysis"""
    print("### PC loadings")
    print("%-17s:\t%10s\t%10s" % ('Feature', 'coef PC1', 'coef PC2'))
    for ft, coef1, coef2 in sorted(zip(features, pca.components_[0,:],
                                       pca.components_[1,:]),
                                   key=lambda x: (abs(x[1]), abs(x[2])),
                                   reverse=True):
        print("%-17s:\t%10.6f\t%10.6f" % (ft, coef1, coef2))


def detailed_pca(alls_normed, features, FA=False, abs_cov=True, make_corr=True,
                 heat_dendro=True, out=None, **heat_kwargs):
    """if FA=True, perform factor analysis instead of PCA"""

    disp = getattr(out, 'output', display)
    disp_html = getattr(out, 'html', display_html)
    show = getattr(out, 'show', plt.show)  # Bad things may happen with unexpected classes.

    Analysis = FactorAnalysis if FA else PCA
    fa = Analysis(n_components=15)
    transformed = fa.fit_transform(alls_normed[features])

    outputs = []
    # ### PC contribution to variance

    print("transformed data dimensions: %s" % (transformed.shape,), file=out)
    print("components dimensions: %s" % (fa.components_.shape,), file=out)

    # Explained variance of each principal component.
    if FA:
        print('WARNING: Computing feature contribution from FactorAnalysis without the noise variance!', file=out)
        #PC_contrib = (fa.components_**2 + fa.noise_variance_).sum(axis=1)
        PC_contrib = (fa.components_**2).sum(axis=1)
        PC_contrib /= PC_contrib.sum()
    else:
        PC_contrib = fa.explained_variance_ratio_

    print("Feature contributions:\n", PC_contrib, file=out)
    # Plot cumulative contribution of PC
    fig, ax = plt.subplots()
    ax.bar(np.arange(PC_contrib.size), PC_contrib.cumsum())
    ax.set_title("Cumulative ratio of variance explained by the Principal Components")
    ax.set_xlabel("Principal Components")
    ax.set_ylabel("Cumulative ratio of variance explained")
    show()
    outputs.append(fig)
    plt.close(fig)  # Not sure of the consequences.

    # ### PC loadings

    # weights of each feature in the PCs

    components = pd.DataFrame(fa.components_.T, index=features,
                              columns=["PC%d" % (i+1) for i in
                                       range(fa.components_.shape[0])])

    # Just so that the features with similar component values are grouped together.
    # TO give more weight in the clustering to the first PC, multiply the PC
    # by their eigenvalue.
    ft_dendro = hclust.dendrogram(
                            hclust.linkage(components*PC_contrib,
                                           'average'),
                            labels=features,
                            count_sort='descending',
                            no_plot=True)
    ordered_ft = ft_dendro['ivl']
    ft_order = ft_dendro['leaves']

    #.sort_values(["PC1", "PC2"])
            #apply(centered_background_gradient, cmap="PRGn", extend=0.15).\
    styled_components = components.loc[ordered_ft].style.\
            background_gradient(cmap="PRGn", low=0.85, high=0.85).\
            set_caption("Principal Components loadings").\
            set_properties(**{'max-width': '60px', 'font-size': '1pt'}).\
            set_table_styles(magnify())
    #print("Rendered_components:", type(styled_components), file=out)
    disp_html(styled_components)
    outputs.append(styled_components)

    # ### Feature covariance
    ft_cov = fa.get_covariance()
    print("Covariance dimensions:", ft_cov.shape, file=out)
    if abs_cov:
        ft_cov = np.abs(ft_cov)
    if make_corr:  # heat_corr
        ft_cov = cov2cor(ft_cov)
    if heat_dendro:
        fig = heatmap_cov(ft_cov, features, make_corr=False, **heat_kwargs)
    else:
        plot_cov(ft_cov[ft_order,:][:,ft_order], ordered_ft, cmap='seismic')
        fig = plt.gcf()

    fig.suptitle('Feature %s%s (%s)' % (
        'absolute ' if abs_cov else '',
        'correlation' if make_corr else 'covariance',
        'Factor Analysis' if FA else 'PCA'))
    show()
    outputs.append(fig)
    plt.close(fig)

    # Plot feature vectors in the PC space
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 20),
                                   subplot_kw={'projection': 'polar'})

    plot_features_radar(components, features, ax=ax0)
    plot_features_radar(components, features, PCs=["PC1", "PC3"], ax=ax1)
    fig.suptitle("Features in Principal Component space (%s)" % ('Factor Analysis' if FA else 'PCA'))
    show()
    outputs.append(fig)
    plt.close(fig)

    fig, (ax0, ax1) = plt.subplots(1, 2, sharey=True)
    scatter_density(transformed[:,1], transformed[:,0], alpha=0.5, ax=ax0)
    ax0.set_xlabel('PC2')
    ax0.set_ylabel('PC1')
    scatter_density(transformed[:,2], transformed[:,0], alpha=0.5, ax=ax1)
    ax1.set_xlabel('PC3')
    ax1.set_ylabel('PC1')
    fig.suptitle('Data plotted in component space (%s)' % ('Factor Analysis' if FA else 'PCA'))
    show()
    outputs.append(fig)
    plt.close(fig)
    return fa, outputs


# Functions for checking colinearity between variables

def display_decorrelate_attempts(var, correlated_var, data, logdata=None):
    nrows = 2
    ncols = 1 if logdata is None else 2
    fig = plt.figure()
    
    ax0 = fig.add_subplot(nrows, ncols, 1)
    ax1 = fig.add_subplot(nrows, ncols, 1+ncols, sharex=ax0)
    if logdata is not None:
        ax2 = fig.add_subplot(nrows, ncols, 2)
        ax3 = fig.add_subplot(nrows, ncols, 2+ncols, sharex=ax2)
        
    scatter_density(correlated_var, var, alpha=0.5, s=9, data=data, ax=ax0)
    
    # Plot the decorrelated variables
    #if var not in data.columns and re.search(r'[+-/*]', var):
    #    var = data[var]
    #    correlated_var = data[correlated_var]

    decor_var = data[var] / data[correlated_var]
    decor_null = (data[var] == 0) & (data[correlated_var] == 0)
    decor_var[decor_null] = 0
    scatter_density(data[correlated_var], decor_var, s=9, alpha=0.5, ax=ax1)
    
    raw_na = np.isnan(data[var]) | np.isnan(data[correlated_var])
    decorr_na = np.isnan(decor_var) | np.isnan(data[correlated_var])

    print('%d NAs in original data; %d in decorrelated.' % (raw_na.sum(), decorr_na.sum()))
    print('Pearson R: %.4f -> %.4f' % (
          stats.pearsonr(data[var][~raw_na],
              data[correlated_var][~raw_na])[0], 
          stats.pearsonr(decor_var[~decorr_na],
              data[correlated_var][~decorr_na])[0]))
    print('Spearman R: %.4f -> %.4f' % (
          stats.spearmanr(data[var][~raw_na],
               data[correlated_var][~raw_na])[0], 
          stats.spearmanr(decor_var[~decorr_na],
               data[correlated_var][~decorr_na])[0]))
    ax0.set_title("Original scale")
    ax0.set_ylabel("Original %s" % var)
    ax1.set_ylabel("Decorrelated %s" % var)
    ax1.set_xlabel(correlated_var)
    
    if logdata is not None:
        # Plot the log of decorrelated variables
        logdecor_var = logdata[var] / logdata[correlated_var]
        scatter_density(correlated_var, var, s=9, alpha=0.5, data=logdata, ax=ax2)
        scatter_density(logdata[correlated_var], logdecor_var, s=9, alpha=0.5, ax=ax3)
        #ax2.set_xlabel("log(%s)" % correlated_var)
        ax3.set_xlabel("log(%s)" % correlated_var)
        ax2.set_title("Log-scale")
        print('Pearson R (log VS log): %.4f -> %.4f' %(
              stats.pearsonr(logdata[var], logdata[correlated_var])[0], 
              stats.pearsonr(logdecor_var, logdata[correlated_var])[0]))


def display_decorrelate(decorr_item, data_raw, data_transformed, data_decorred,
                        color_var=None, decorred_corrvar=None, alpha=0.4):
    if color_var is None:
        def scatter(corrvar, yvar, data, ax):
            return scatter_density(corrvar, yvar, data=data, ax=ax, alpha=alpha)
    else:
        #cnorm = mpl.colors.Normalize(data_raw[color_var].min(), data_raw[color_var].max())
        # No, different data have different transforms.
        color_values = data_raw[color_var]
        ncolors = len(data_raw[color_var].unique())
        if np.issubdtype(data_raw[color_var].dtype, np.integer) and ncolors<=10:
            if ncolors == 2:
                cmap = 'PRGn_r'
            else:
                cmap = plt.get_cmap('tab10', len(data_raw[color_var].unique()))
        else:
            cmap = 'viridis'
        def scatter(corrvar, yvar, data, ax):
            try:
                data_s = data.sort_values(color_var)  # for plotting order same as value
            except TypeError as err:
                err.args = (err.args[0] + " in %r" % color_var,)
                raise
            return ax.scatter(corrvar, yvar, data=data_s, c=color_var, alpha=alpha,
                              label=None, cmap=cmap)
    k, (var, corrvar) = decorr_item
    fig, axes = plt.subplots(ncols=3)
    for ax, data, desc in zip(axes, (data_raw, data_transformed, data_decorred),
                              ('raw', 'transformed', 'decorred')):
        if desc=='decorred':
            #FIXME: This will raise errors if no decorrelation was done!
            # Alt solution: provide k == yvar
            if yvar not in data:
                yvar = 'R'+var if isinstance(k, int) else k
            else:
                logger.warning('%r was not decorred!', yvar)

            if decorred_corrvar is None:
                # The correlated variable was left unchanged.
                data = data.assign(**{corrvar: data_transformed[corrvar]})
            else:
                # The correlated variable was also decorred against another one!
                # It's the last data, so we can update corrvar:
                corrvar = decorred_corrvar

        else:
            yvar = var
        points = scatter(corrvar, yvar, data, ax)
        ax.set_title(desc)
        ax.set_xlabel(corrvar)
        ax.set_ylabel(yvar)
        if color_var is not None:
            fig.colorbar(points, ax=ax, orientation='horizontal', fraction=0.05).set_label(color_var) #aspect=30

        if (data[[yvar, corrvar]].isna().any(axis=None) or
                np.isinf(data[[yvar, corrvar]]).any(axis=None)):
            n_NaNs = data[[yvar, corrvar]].isna().sum().values
            n_Infs = np.isinf(data[[yvar, corrvar]]).sum().values
            logger.error('NaN: %s, Inf: %s in %s [%s, %s]', n_NaNs, n_Infs, desc, yvar, corrvar)
            data = data.loc[np.isfinite(data[[yvar, corrvar]]).all(axis=1)]
            if not data.shape[0]:
                logger.error('ZERO valid rows in %s [%s, %s]!', desc, yvar, corrvar)
                continue

        fit = sm.OLS(data[yvar], sm.add_constant(data[corrvar])).fit()
        if desc=='raw':
            # xlims are set to large
            xmin, xmax = data[corrvar].min(), data[corrvar].max()
            xrng = xmax - xmin
            xmin, xmax = xmin - xrng*0.05, xmax + xrng*0.05
            ax.set_xlim(xmin, xmax)
            ymin, ymax = data[var].min(), data[var].max()
            yrng = ymax - ymin
            ax.set_ylim(ymin - yrng*0.05, ymax + yrng*0.05)
        else:
            xmin, xmax = ax.get_xbound()
        x = np.linspace(xmin, xmax, 2)
        a, b = fit.params
        r2, pval = fit.rsquared, fit.f_pvalue
        ax.plot(x, a + b*x, '--', alpha=alpha,
                label='a=%g b=%g\nR²=%.3f P-v=%.3g' % (a, b, r2, pval))
        ax.legend()
    #fig.tight_layout()  # tight layout often goes havoc with colorbars.
    return fig


def display_regress_scatter(xvar, yvar, data=None, ax=None):
    scatter_density(xvar, yvar, alpha=0.3, s=9, data=data, ax=ax)
    if ax is None: ax = plt.gca()
    xlim = np.array(ax.get_xlim())
    ylim = ax.get_ylim()

    #if var not in data.columns and re.search(r'[+-/*]', var):
    #    var = data[var]
    #    correlated_var = data[correlated_var]
    if data is None:
        xdata = xvar
        ydata = yvar
        try:
            xvar = xdata.name
        except AttributeError:
            xvar = 'x'
            xdata = pd.Series(xdata, name='x')
    else:
        xdata = data[xvar]
        ydata = data[yvar]

    finite = np.isfinite(xdata) & np.isfinite(ydata)
    xdata = xdata[finite]
    ydata = ydata[finite]

    fit = sm.OLS(ydata, sm.add_constant(xdata)).fit()
    a, b = fit.params[['const', xvar]]
    #ax.annotate('', xright, ybottom, ''
    ax.plot(xlim, a + b*xlim, '-', label='y = %g + %g x\nR²=%.3f  F p-val=%.3g' % (a,b,fit.rsquared, fit.f_pvalue))
    ax.annotate(('%d non finite\n' % (~finite).sum()) +
        ('Pearson R: %.4f\n' % stats.pearsonr(xdata, ydata)[0]) +
        ('Spearman R: %.4f\n' % stats.spearmanr(xdata, ydata)[0]),
        (xlim[1], ylim[1] if b<0 else ylim[0]),
        ha='right', va=('top' if b<0 else 'bottom'))
    ax.legend()
    return ax



def bootreg(y, x, data, n=100, model=sm.OLS, fitmethod="fit", add_const=True, **fitkw):
    #data = fit.model.data
    if isinstance(x, str):
        x = [x]
    if add_const:
        data = sm.add_constant(data)  # Check 'const' is in x
        if 'const' not in x:
            x.append('const')
    coeffs = np.empty((len(x), n))
    for i in range(n):
        newdata = data.sample(frac=1, replace=True)
        reg = model(newdata[y], newdata[x], hasconst=True)
        fit = getattr(reg, fitmethod)(**fitkw)
        # WARNING: because of sampling with replacement, X might not be invertible
        #print(
        coeffs[:,i] = fit.params.values
    return coeffs

def bootreg_summary(coeffs, alpha=0.05):
    return (coeffs.mean(axis=1),
            coeffs.std(ddof=0, axis=1),
            np.percentile(coeffs, [alpha/2., 1-alpha/2.], axis=1))


def bootreg_fromfit(fit, n=100):
    y = fit.model.endog_names
    x = fit.model.exog_names
    model = type(fit.model)
    raise NotImplementedError
    #if isinstance(fit, statsmodels.base.elastic_net.RegularizedResultsWrapper):
    #    fitmethod = 'fit_regularized'
    #else:
    #    fitmethod = 'fit'
    #return bootreg(y, x,



def pairwise_regress_stats(data, features):
    """Return R², intercept, slope, Pval."""
    N = len(features)
    # Attributes and indices to get the statistics:
    stats = {'R2': ('rsquared',None),
             'const': ('params', 0),
             'slope': ('params', 1),
             'Pval': ('f_pvalue', None)}
    mats = {stat: np.full((N, N), np.NaN) for stat in stats}
    for j, fty in enumerate(features[:-1]):
        for i, ftx in enumerate(features[j+1:], start=j+1):
            fit = sm.OLS(data[fty], sm.add_constant(data[[ftx]])).fit()
            for stat, (attr, idx) in stats.items():
                value = getattr(fit, attr)
                if idx is not None:
                    value = value[idx]
                mats[stat][j, i] = value  # Upper triangular.
    # TODO: add the 1 against 2 combinations
    return mats


def leave1out_eval(features, func, dtype=np.float):
    """Example: func=lambda x: multicol_test(df[x])"""
    out = pd.Series(0, index=features, dtype=dtype)
    for i, ft in enumerate(features):
        out[features[i]] = func(features[:i] + features[(i+1):])

    return out


def display_leave1out_eval(features, func):
    return leave1out_eval(features, func).to_frame.style.bar()


symbol_to_operator = {'<=': '__le__',
                      '>=': '__ge__',
                      '==': '__eq__',
                      '!=': '__ne__',
                      '<': '__lt__',
                      '>': '__gt__'}
stop_reg = re.compile(r'(' + r'|'.join(symbol_to_operator) + r')\s*(.*)')

def loop_leave1out(features, func, criterion='min', nloops=None, stop_criterion=None,
                   protected=None, start_value=None,
                   na_equiv='nan', format_df='matplotlib',
                   widget=None, out=None):
    """
    `criterion`: how to select what to *drop* ("min"/"max").
    `stop_criterion`: string with operator and number: example: ">20".
    
    NOTE: NaN values returned by `func` never fulfill the stop_criterion.
    na_equiv = "nan" -> NA values propagate in min and max.
    na_equiv = "Inf" -> NA values propagate in max only.
    na_equiv = "-Inf" -> NA values propagate in min only.
    na_equiv = "ignore" -> NA values are skipped.

    EXAMPLE:
        func      = lambda x: multicol_test(df[x]),
        criterion = "min",
        na_equiv  = "max"  # because multicol_test returns NaN when super high colinearity!
    """
    #if out is not None:
        # Replace the functions: `print`, `display`, `plt.show`
        #FIXME: Implementation could probably be nicer with something like Mock.patch.
    disp = getattr(out, 'output', display)
    disp_html = getattr(out, 'html', display_html)
    show = getattr(out, 'show', plt.show)  # Bad things may happen with unexpected classes.

    protected = [] if protected is None else list(protected)

    if nloops is None:
        nloops = len(set(features).difference(protected))

    outputs = []
    if format_df is None:
        def format_df(series):
            pass
    elif format_df in ('mpl', 'matplotlib'):
        rcwidth, rcheight = mpl.rcParams['figure.figsize']
        bw = max(0.15, rcwidth/len(features))
        bh = min(rcheight, bw*5)
        def format_df(series):
            ax, = matplotlib_stylebar(series, horizontal=False, ticks=True)
            ax.figure.set_size_inches(bw*series.shape[0], bh)
            show()
            outputs.append(ax.figure)
            plt.close(ax.figure)
    elif format_df in ('pd', 'pandas'):
        def format_df(df):
            styled = series.to_frame().style.bar()
            disp_html(styled)
            outputs.append(styled)
    elif not callable(format_df):
        raise ValueError('`show` must be a value in [None, "pandas", "matplotlib"] or a callable')

    if widget is None:
        widget = lambda iterator, *a, **kw: iterator
    
    na_equiv = na_equiv.lower()
    if criterion == 'min':
        select_drop = np.argmin if na_equiv in ('nan', '-inf') else np.nanargmin
    elif criterion == 'max':
        select_drop = np.argmax if na_policy in ('nan', 'inf') else np.nanargmax
    else:
        raise ValueError('Unknown criterion %r (min/max)' % criterion)
    if stop_criterion is None:
        def is_stopped(values, i): return False
    else:
        # Old behavior: raise DeprecationWarning
        if isinstance(stop_criterion, (int, float)):
            if criterion == 'min':
                def is_stopped(values, i): return values[i] < stop_criterion
            elif criterion == 'max':
                def is_stopped(values, i): return values[i] > stop_criterion
        else:
            stop_comp, stop_val = stop_reg.match(stop_criterion).groups()
            stop_comp = symbol_to_operator[stop_comp]
            try:
                stop_val = int(stop_val)
            except ValueError:
                stop_val = float(stop_val)
            
        def is_stopped(values, i): return getattr(values[i], stop_comp)(stop_val)

    dropped_features = []
    if is_stopped([start_value if start_value is not None else func(features)], 0):
        return dropped_features, outputs

    next_features = [ft for ft in features]
    for k in widget(range(nloops)):
        left1out_k = leave1out_eval(next_features, func)
        if na_equiv == 'nan' and np.isnan(left1out_k[protected].values).any():
            logger.warning('NaN values in func evaluations (protected features).')
        format_df(left1out_k)
        # Also check protected features
        i = select_drop(left1out_k.values)
        if next_features[i] in protected:
            logger.info('Would have dropped protected feature %r',
                        left1out_k)
            # "Protect" from criterion by setting a "good" value (min if "max" and conversely)
            left1out_k[protected] = left1out_k.min() if criterion=='max' else left1out_k.max()
            i = select_drop(left1out_k.values)

        if na_equiv == 'nan' and np.isnan(left1out_k.drop(protected).values).any():
            logger.warning('NaN values in func evaluations of (non-protected features).')
        stopped = is_stopped(left1out_k.values, i)
        
        if next_features[i] in protected:
            # It means that there are no un-protected features to drop. Must exit.
            break
        print('%d. DROP %s' % (k, next_features[i]), file=out)
        dropped_features.append(next_features[i])
        
        if stopped:
            break
        next_features = next_features[:i] + next_features[(i+1):]

    if not stopped:
        logger.warning('could not reach the criterion %s', stop_criterion)
    return dropped_features, outputs


# Functions for linear models

def lm_summary(lm, features, response, data):
    """Display the summary statistics of the sklearn multiple linear regression."""
    print("R^2       =", lm.score(data[features], data[response]))
    print("Intercept =", lm.intercept_)
    print("\nSlopes\n------")

    features_by_coef = sorted(zip(features, lm.coef_), key=lambda x: np.abs(x[1]), reverse=True)

    for ft, coef in features_by_coef:
        print("%-17s: %10.6f" % (ft, coef))


PVAL_REGEX = re.compile(r'p>\|z\||p-?val(ue)?', re.I)


def sm_pretty_slopes(olsfit, join=None, renames=None, bars=['coef'], bars_mid=None):
    """Nicer look for a StatsModels OLS fit (sorted features by slope)."""
    try:
        summary = olsfit.summary()
    except NotImplementedError:
        summary = None
    if summary is not None:
        r_coefs = pd.read_csv(StringIO(summary.tables[1].as_csv()),
                              sep=r'\s*,\s*', index_col=0, engine='python')
    else:
        # For example a Lasso fit (OLS().fit_regularized())
        r_coefs = olsfit.params.to_frame('coef')

    # Sort
    coef_order = r_coefs.coef.abs().sort_values(ascending=False).index
    if 'const' in coef_order:
        coef_order = ['const'] + coef_order.difference(('const',), sort=False).tolist()

    r_coefs = r_coefs.loc[coef_order]

    if join is not None:
        # Add additional info about the params.
        r_coefs = r_coefs.join(join,sort=False)

    renames = {} if renames is None else renames

    param_names = [renames.get(n, n) for n in olsfit.model.exog_names
                   if n not in ('Intercept', 'const')]
    r_coefs_styled = r_coefs.rename(renames)\
                     .style.bar(
                        subset=pd.IndexSlice[param_names, bars],
                        axis=0,
                        align='zero')
    if bars_mid:
        r_coefs_styled = r_coefs_styled.bar(subset=pd.IndexSlice[param_names, bars_mid],
                                            axis=0, align='mid')
    pval_columns = [col for col in r_coefs.columns if PVAL_REGEX.search(col)]
    if pval_columns:
        r_coefs_styled = r_coefs_styled.applymap(
                         lambda v: ('background: khaki' if v<0.01 else
                                    'background: lemonchiffon; alpha: 0.5' if v<=0.05 else ''),
                                 subset=pval_columns)
    return r_coefs_styled


def sm_pretty_summary(fit, join=None, renames=None, bars=['coef'], bars_mid=None,
                      out=None):
    try:
        summary = fit.summary()
    except NotImplementedError:
        summary = None
    if summary is not None:
        try:
            out.html(summary.tables[0])
            for table in summary.tables[2:]:
                out.html(table)
        except AttributeError:
            display_html(summary.tables[0])
            for table in summary.tables[2:]:
                display_html(table)
    else:
        print('Warning: R² and adjusted R² not provided for this result.', file=out)
        print('R² =', r_squared(fit.model.endog, fit.fittedvalues),
              '; Adj. R² =', adj_r_squared(fit.model.endog,
                                         fit.fittedvalues,
                                         len(fit.params)), file=out)
    pretty_slopes = sm_pretty_slopes(fit, join, renames, bars, bars_mid)
    #display_html(pretty_slopes)
    return pretty_slopes

## How to use sm_pretty_summary with HtmlReport:
#with HtmlReport('outfile.html') as hr:
#    with contextlib.capture_stdout(hr):
#        pretty_slopes = sm_pretty_summary()
#    hr.html(pretty_slopes)
