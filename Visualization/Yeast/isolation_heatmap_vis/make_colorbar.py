import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
mpl.rcParams["axes.labelcolor"] = "#000000"
mpl.rcParams["axes.linewidth"]  = 1.0 

def generate_cmap(colors):
    color_list = []
    values = range(len(colors))
    vmax   = int(np.max(values))
    for v, c in enumerate(colors):
        color_list.append( (v*1.0/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


fig = plt.figure(figsize=(4,4))

# define the data
x = np.random.rand(10)
y = np.random.rand(10)
tag = np.random.randint(0,10,10)
tag[10:12] = 0 # make sure there are some 0 values to showup as grey

#cmap     = plt.cm.Blues
#cmaplist = [cmap(i) for i in range(0,256)][0:]
#cmap     = generate_cmap(["#FFFFE2",cmaplist[-20]])
#cmap = plt.cm.YlGn
#cmap = [cmap(i) for i in range(0,256)][0:]
#cmap = plt.cm.PuBu
#cmap = [cmap(i) for i in range(0,256)][0:]
cmap = plt.cm.Wistia
#cmap = [cmap(i) for i in range(0,256)][0:]


ax2 = fig.add_axes([0.95, 0.1, 0.06, 0.36])
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap)
cb.ax.tick_params(length=3,which="minor")
cb.set_ticks([0,0.5,1.0])
cb.set_ticklabels([])
ax2.tick_params(labelsize=20)
fig.patch.set_alpha(0.0)
fig.savefig("hoge.pdf",bbox_inches="tight")
