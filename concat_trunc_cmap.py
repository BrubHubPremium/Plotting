import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

interval = np.linspace(.5,.75,128)
colors11 = plt.cm.seismic(interval)
colors2 = plt.cm.seismic(np.linspace(0, 1, 128))
colors1 = plt.cm.autumn(np.linspace(0, 1, 128))
colors = np.vstack((colors11, colors1))
gangstacmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
