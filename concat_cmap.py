import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
# Positive | Negative (colors1 | colors2)
colors1 = plt.cm.Greens(np.linspace(0, 1, 128))
colors2 = plt.cm.Oranges(np.linspace(1, 0, 128))
colors = np.vstack((colors2, colors1))
mymap = mcolors.LinearSegmentedColormap.from_list(‘my_colormap’, colors)
