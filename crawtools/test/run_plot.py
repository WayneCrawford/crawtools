import sys
from crawtools.bathy_map import BathyMap
from pylab import cm

gridx, gridy = 0.05, 0.05
bounds = [-32.4, -32.2, 37.2, 37.4]
# bounds = [-33, -31.8, 36.8, 38]
bathy_file = '/Users/crawford/_Work/Figures_Etc/2_Maps/LuckyStrike/grd/LuckySISMOMAR_40mtr.grd'
plot_filename = 'plot_bathy_40m.png'
plot_title = 'East pressure monitoring site'

b = BathyMap(bounds, bathy_file, gridx, gridy)
b.plot_image(cmap=cm.Spectral_r, scale=1, azdeg=-90.0, altdeg=45.0)
b.plot_contours(500, linewidth=0.5)
b.plot(-32.2814, 37.2926, color='red') # N37d17.559 / W32d16.885
b.text(-32.2814 + 0.01, 37.2926, 'JPPW', ha='left', va='center') # N37d17.559 / W32d16.885
b.plot(-32.2475, 37.2831) #N37d16.986 / W32d14.853
b.text(-32.2475 + 0.01, 37.2831, 'JPPE', ha='left') #N37d16.986 / W32d14.853
b.save_map(plot_filename, plot_title)
#b.savefig(plot_filename, dpi=300)