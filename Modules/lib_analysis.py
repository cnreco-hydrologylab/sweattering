"""
Library: lib_analysis.py
Author(s): Martina Natali (martinanatali@cnr.it)
Date: 2024-02-19
Version: 1.0.0
"""
# Base
import math

# Graphics
import matplotlib.pyplot as plt
import numpy as np

# Analysis
from scipy import special as sp
from scipy import stats
from scipy.optimize import curve_fit


#############################################################################
# Fitting functions
#############################################################################

def gauss(x, sigmag, mu, A):
    """Not-normalized, shifted gaussian distribution."""
        
    pdf = (1/(sigmag*np.sqrt(2*math.pi)))*np.exp(-(x-mu)**2/(2*sigmag**2))
    return A*pdf


def skew_gauss(x, sigmag, mu, alpha, a):
    
    normpdf = (1/(sigmag*np.sqrt(2*math.pi)))*np.exp(-(np.power((x-mu),2)/(2*np.power(sigmag,2))))
    normcdf = (0.5*(1+sp.erf((alpha*((x-mu)/sigmag))/(np.sqrt(2)))))
    
    return 2*a*normpdf*normcdf


#############################################################################
# Histogram with fit
#############################################################################

def hist_gauss_fit(data, nbins, hist_kwargs, fitline_kwargs,
                   title, density=False,
                   opt_save=False, dir_name='', opt_name='hist_fit',
                   opt_fit=True,
                   func='optim', thr_asymm=2, fit_method='trf',
                  ):
    """Histogram with automatic gaussian fit.
    
    Arguments
    ---------
    - func: object, default gauss
        WARNING: skew_gauss not supported yet
        
    """
    
    from scipy import integrate

    fontdict = {
        'size':13
    }
    
    mean = np.mean(data); rang=np.ptp(data);
    x = np.linspace(mean-2*rang, mean+2*rang, 500)
    counts, bins, pads = plt.hist(data, bins=nbins, density=density, **hist_kwargs)
    
    thr_asymm = thr_asymm
    if func=='optim':
        if abs(stats.skew(data))<thr_asymm: func_fit=gauss
        else: func_fit=skew_gauss
    else: func_fit=func
        
    integ = sum(np.diff(bins)*counts)/sum(counts) if density==False else 1
    
    if func_fit==gauss:
        fit_label='Gaussian'
        kwargs = {
        'bounds': [ 
            # sigmag, mu, a
            [0,min(data),0],
            [np.std(data)*2,max(data),integ*2] ]
        }
    elif func_fit==skew_gauss:
        fit_label='Skew normal'
        kwargs = {
            'p0' : (np.std(data), np.mean(data),
                    stats.skew(data), integ),
            'bounds' : [
                # sigmag, mu, alpha, a
                [0, min(data),  min(0,stats.skew(data)*2), min(0,integ*2),],
                [np.std(data)*2, max(data), max(0,stats.skew(data)*2), max(0,integ*2),]
            ]
        }
    else: raise ValueError(f'Func {func} is not a valid option.')

    
    if opt_fit==True:
        options = {
            # 'loss' : 'cauchy',
        }

        popt, pcov = curve_fit(func_fit, bins[:-1], counts,
                               method=fit_method, **kwargs, **options,
                               maxfev=100000, #full_output=False
                               )
        fit = func_fit(x, *popt)
        plt.plot(x, fit, label=fit_label+' fit', **fitline_kwargs)
        mode = x[np.argmax(fit)]
        
        def func_fit_popt(x): return func_fit(x, *popt)
        
        def quantiles(x):
            integral = integrate.quad(func_fit_popt, -np.inf, np.inf)[0]
            pdf = lambda x: func_fit(x, *popt)/integral
            cdf = lambda x: integrate.quad(pdf, -np.inf, x)[0]
            y = np.array([]); q1 = False; median = False; q3 = False
            
            for xi in x:
                yi = cdf(xi)
                if (not q1)and(yi>=1/4): q1=xi
                if (not median)and(yi>=1/2): median=xi
                if (not q3)and(yi>=3/4): q3=xi
                y = np.append(y, yi)

            return [q1, median, q3]
        
        q1, q2, q3 = quantiles(x)

    else:
        fit=[]; popt = []; pcov = []
        data_sorted = np.sort(data)
        q1 = data_sorted[int(len(data)*1/4)]
        q2 = data_sorted[int(len(data)*2/4)]
        q3 = data_sorted[int(len(data)*3/4)]
        mode = bins[np.argmax(counts)]
        
        
    maxy = np.max([*fit, *counts])
    plt.vlines(q1, 0, maxy, color='k', linestyle='--', alpha=.7 , label='Q1')
    plt.vlines(q2, 0, maxy, color='k', linestyle='-', alpha=1 ,  label='Q2')
    plt.vlines(q3, 0, maxy, color='k', linestyle='--', alpha=.7 , label='Q3')

    ylabel = 'Density' if density else 'Counts';
    plt.ylabel(ylabel, fontdict=fontdict)
    plt.xlabel(title, fontdict=fontdict);
    plt.legend(loc='best', framealpha=.5);
    mean=np.mean(data); delta=15*(q3-q1); rang=np.ptp(data)
    plt.xlim(q2-min([delta, rang]), q2+min([delta, rang]));
    
    plt.tick_params(axis='both', which='both', right=False, top = False)
    plt.tick_params(axis='both', which='both', labelsize=14)
    
    if opt_save: plt.savefig(dir_name+opt_name+'.png', dpi=300)
    
    
    return [counts, bins, pads, popt, pcov, q1, q2, q3, mode]


#############################################################################
# Analysis on parameters out of multiple runs
#############################################################################

def parameters_analysis(params, PAR_dict, PAR, PARn_str, bounds, timestr,
                        automate, nbins=10, opt_fit=True, root='/Plot'):
    matrix = np.array(
        [
            np.array(
                [params[i].take(j) for i in range(len(params))])
            for j in range(len(PAR))
        ]
    )

    PARn = []
    PARn_dev = []
    rowLabels = []

    if automate:
        opt_save = True
    else:
        opt_save = True if input(
            'Save histograms of params? [y/n]') == 'y' else False
    opt_show = True

    for label in PAR_dict:
        if label in PARn_str:
            i = PARn_str.index(label)
            data = matrix[i]

            if len(data) == 1:
                PARn = [x[0] for x in matrix]
                PARn_dev = [0 for x in matrix]
            else:
                hist_kwargs = {'alpha': .5, }
                fitline_kwargs = {'linestyle': '-', }
                counts, bins, pads, popt, pcov, q1, q2, q3, mode = hist_gauss_fit(
                    data, nbins=nbins, hist_kwargs=hist_kwargs,
                    fitline_kwargs=fitline_kwargs,
                    title=f'{label} {PAR_dict[label][0]}', density=True,
                    opt_save=opt_save, dir_name=root,
                    opt_name=f'{timestr}_hist_{i}',
                    opt_fit=opt_fit,
                    func=gauss, fit_method='dogbox')
                mean = np.mean(data);
                rang = np.ptp(data)
                plt.xlim(mean - rang, mean + rang);

                PARn.append(q2)
                PARn_dev.append((q3 - q1) / 2)
                plt.show()

            if opt_show:
                plt.show()
            else:
                plt.close()

            # PAR_dict[label].pop(1)
            PAR_dict[label].append('cal');
            PAR_dict[label].append(
                [int(round(bounds[0][i], 2)), int(round(bounds[1][i], 2))]);
            PAR_dict[label].append(round(PARn[i], 3));
            PAR_dict[label].append(round(PARn_dev[i], 3));
            PAR_dict[label].append(list(matrix[i]));

            if not 'Guess' in rowLabels: rowLabels.append('Guess')
            if not 'Cal/fix' in rowLabels: rowLabels.append('Cal/fix')
            if not 'Bounds' in rowLabels: rowLabels.append('Bounds')
            if not 'Median' in rowLabels: rowLabels.append('Median')
            if not 'Err' in rowLabels: rowLabels.append('Err')

        else:
            PAR_dict[label].append('fix')
            PAR_dict[label].append(['/', '/'])
            PAR_dict[label].append(PAR_dict[label][1])
            PAR_dict[label].append('/')
            PAR_dict[label].pop(1)

    # print(PAR_dict, '\n', rowLabels)
    return PAR_dict, rowLabels