# Create plots for ICC results:
import os,sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
def _bivariate_kdeplot(xx, yy, z, filled, fill_lowest,
                       kernel, bw, gridsize, cut, clip,
                       axlabel, cbar, cbar_ax, cbar_kws, ax, **kwargs):
    from seaborn.palettes import color_palette, light_palette, dark_palette, blend_palette
    """Plot a joint KDE estimate as a bivariate contour plot."""
    # Determine the clipping
    if clip is None:
        clip = [(-np.inf, np.inf), (-np.inf, np.inf)]
    elif np.ndim(clip) == 1:
        clip = [clip, clip]
    # Plot the contours
    n_levels = kwargs.pop("n_levels", 10)
    scout, = ax.plot([], [])
    default_color = scout.get_color()
    scout.remove()
    color = kwargs.pop("color", default_color)
#     cmap = kwargs.pop("cmap", None)
    cmap = cbar_kws['cmap']
#     if cmap is None:
#         if filled:
#             cmap = light_palette(color, as_cmap=True)
#         else:
#             cmap = dark_palette(color, as_cmap=True)
#     if isinstance(cmap, string_types):
#         if cmap.endswith("_d"):
#             pal = ["#333333"]
#             pal.extend(color_palette(cmap.replace("_d", "_r"), 2))
#             cmap = blend_palette(pal, as_cmap=True)
#         else:
#             cmap = mpl.cm.get_cmap(cmap)
    label = kwargs.pop("label", None)
    kwargs["cmap"] = cmap
    contour_func = ax.contourf if filled else ax.contour
    cset = contour_func(xx, yy, z, n_levels, **kwargs)
    if filled and not fill_lowest:
        cset.collections[0].set_alpha(0)
    kwargs["n_levels"] = n_levels
    if cbar:
        cbar_kws = {} if cbar_kws is None else cbar_kws
        ax.figure.colorbar(cset, cbar_ax, ax, **cbar_kws)
    # Label the axes
    if hasattr(xx, "name") and axlabel:
        ax.set_xlabel(x.name)
    if hasattr(yy, "name") and axlabel:
        ax.set_ylabel(y.name)
    if label is not None:
        legend_color = cmap(.95) if color is None else color
        if filled:
            ax.fill_between([], [], color=legend_color, label=label)
        else:
            ax.plot([], [], color=legend_color, label=label)
    return ax
# ex: flank 5sessions
# task = sys.argv[1]
# subset = sys.argv[2]
ratios = {}
for task in ['rest','flank','inscape','movie']:
    subset = '2subset'

    csvdir = '/data3/cdb/jcho/ssihyp/analysis/reliability/icc_csv/%s_%s' % (task, subset)
    lhparcels = pd.read_csv('/data3/cdb/anna/resources/parcel/parcel2network/181Yeo7matchlh.csv').values[1:,2]
    rhparcels = pd.read_csv('/data3/cdb/anna/resources/parcel/parcel2network/181Yeo7matchrh.csv').values[1:,2]
    allparcels = np.r_[lhparcels,rhparcels]
    results = {}
    mats = {}

    for metric in ['icc', 'within', 'between', 'vartotal']:
        results[metric] = pd.read_csv('%s/%s_%s_vector.csv' % (csvdir, task, metric)).values[:,1]
        mats[metric] = np.zeros([360,360])

        mats[metric][np.triu_indices(360,k=1)] = results[metric]
        mats[metric] = mats[metric].T
        mats[metric][np.triu_indices(360,k=1)] = results[metric]
        a = mats[metric][np.argsort(allparcels),:]
        a = a[:,np.argsort(allparcels)]
        if metric == 'icc':
            plt.figure(figsize=(10,10));plt.imshow(a,vmin=0,vmax=1);plt.colorbar();plt.savefig('%s/%s_%s_mat.png' % (csvdir, task, metric))
            plt.figure(figsize=(10,10));sns.distplot(mats['icc'][mats['icc']>0.001],hist=False);plt.xlim([0,1]);plt.savefig('%s/%s_%s_hist.png' % (csvdir, task, metric))
        if metric == 'within':
            plt.figure(figsize=(10,10));plt.imshow(a,vmin=0,vmax=0.04);plt.colorbar();plt.savefig('%s/%s_%s_mat.png' % (csvdir, task, metric))
        if metric == 'between':
            plt.figure(figsize=(10,10));plt.imshow(a,vmin=0,vmax=0.1);plt.colorbar();plt.savefig('%s/%s_%s_mat.png' % (csvdir, task, metric))
        if metric == 'vartotal':
            plt.figure(figsize=(10,10));plt.imshow(a,vmin=0,vmax=0.1);plt.colorbar();plt.savefig('%s/%s_%s_mat.png' % (csvdir, task, metric))

    #(6) make a mask parcels x parcels to keep track of the failing connections on LMM
    mask=1*(mats['icc']<0.001)
    
    # Make kde within/between ratio plots:
    if task == 'rest':
        colork = 'Reds'
    elif task == 'flank':
        colork = 'Blues'
    elif task == 'inscape':
        colork = 'Oranges'
    else:
        colork = 'Greens'
    np.save('%s/%s_failed_mask.npy' % (csvdir, task),mask)
    between_ratio = mats['between']/mats['vartotal']
    between_ratio[between_ratio>1] = np.nan
    between_utri = np.asarray(between_ratio[np.triu_indices(len(between_ratio),1)])
    bmask = np.where(~np.isnan(between_utri)==True)[0]
    within_ratio = mats['within']/mats['vartotal']
    within_ratio[within_ratio>1] = np.nan
    within_utri = np.asarray(within_ratio[np.triu_indices(len(within_ratio),1)])
    wmask = np.where(~np.isnan(within_utri)==True)[0]
    totmask = np.intersect1d(bmask,wmask)
    ratios[task] = {'within':within_utri,'between':between_utri,'totmask':totmask}
#     plt.figure(figsize=(10,10))
#     plt.xlim([0,1])
#     plt.ylim([0,1])
#     plt.plot([1,0],[1,0],color='black',alpha=0.3)
#     for iccline in [0.2,0.4,0.6,0.8]:
#         plt.plot([1,0],[iccline,0],color='black',alpha=0.3)
#         plt.plot([iccline,0],[1,0],color='black',alpha=0.3)
#     sns.kdeplot(within_utri[totmask],between_utri[totmask],shade=True,cmap=colork)
#     plt.ylabel('Inter-variation ratio')
#     plt.xlabel('Intra-variation ratio')
#     plt.savefig('%s/%s_ratio_kde.png' % (csvdir, task))


from scipy.stats import gaussian_kde
from seaborn.distributions import _scipy_bivariate_kde
tasks = ['flank','inscape','movie','rest']
taskcolors = {'flank':'Blues','inscape':'Oranges','movie':'Greens','rest':'Reds'}

for cond1 in tasks:
    for cond2 in [i for i in tasks if i != cond1]:
        
        t1color = taskcolors[cond1]
        t2color = taskcolors[cond2]
        
        print(cond1,cond2,t1color,t2color)
        # Get data and create masks so conditions use same edges
        cond1w = ratios[cond1]['within']
        cond1b = ratios[cond1]['between']
        wmask1 = np.where(~np.isnan(cond1w)==True)[0]
        bmask1 = np.where(~np.isnan(cond1b)==True)[0]
        cond1totmask = np.intersect1d(bmask1,wmask1)

        cond2w = ratios[cond2]['within']
        cond2b = ratios[cond2]['between']
        wmask2 = np.where(~np.isnan(cond2w)==True)[0]
        bmask2 = np.where(~np.isnan(cond2b)==True)[0]
        cond2totmask = np.intersect1d(bmask2,wmask2)

        # Mask b/w and w/in values for each condition
        bothmask = np.intersect1d(cond1totmask,cond2totmask)
        cond1w = cond1w[bothmask]
        cond1b = cond1b[bothmask]
        cond2w = cond2w[bothmask]
        cond2b = cond2b[bothmask]

        # Create 2dkde density in the same way as seaborn:
        bw='scott'
        gridsize=100
        cut=3
        clip = [(-np.inf, np.inf), (-np.inf, np.inf)]
        shade=True
        filled=True
        fill_lowest=False

        xx1, yy1, z1 = _scipy_bivariate_kde(cond1w, cond1b, bw, gridsize, cut, clip)
        xx2, yy2, z2 = _scipy_bivariate_kde(cond2w, cond2b, bw, gridsize, cut, clip)
        xx = xx1-xx2
        yy = yy1-yy2
        z = z1-z2

        shade=True
        vertical=False
        kernel="gau",
        bw="scott"
        gridsize=100
        cut=3
        clip=None
        legend=True
        cumulative=False
        shade_lowest=True
        cbar=True
        cbar_ax=None
        cbar_kws={'cmap':mymap}
        plt.figure(figsize=(10,10))
        ax=plt.gca()
        ax.axes.set_xlim([0,1])
        ax.axes.set_ylim([0,1])

        

        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        # Positive
        exec('colors1 = plt.cm.%s(np.linspace(0,1,128))' % t1color)
        # Negative
        exec('colors2 = plt.cm.%s(np.linspace(0,1,128))' % t2color)
        colors = np.vstack((colors2, colors1))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        ax = _bivariate_kdeplot(xx1, yy2, z, shade, shade_lowest, kernel, bw, gridsize, cut, clip, legend, cbar, cbar_ax, cbar_kws, ax, vmin = -max(np.max(z),abs(np.min(z))), vmax = max(np.max(z),abs(np.min(z))))
        ax.plot([1,0],[1,0],color='black',alpha=0.3)
        for iccline in [0.2,0.4,0.6,0.8]:
            ax.plot([1,0],[iccline,0],color='black',alpha=0.3)
            ax.plot([iccline,0],[1,0],color='black',alpha=0.3)
# plt.savefig('/data3/cdb/jcho/ssihyp/analysis/reliability/%s_%s/%s2%s_ratio_kde.png' % ())
