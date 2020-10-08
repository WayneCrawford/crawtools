from crawtools.bathy_map import BathyMap

gridx, gridy = 0.05, 0.05
bounds = [-32.4, -32.2, 37.2, 37.4]
bathy_file = '/Users/crawford/_Work/Figures_Etc/2_Maps/LuckyStrike/grd/LuckySISMOMAR_40mtr.grd'

b = BathyMap(bounds, bathy_file, gridx, gridy)
b.plot_image()
b.plot_contours(500, linewidth=1)
b.plot_contours(100, linewidth=0.5)
# b.plot_coastlines()
b.plot_station(-32.27, 37.25, 'LSVVVVV')
b.plot_station(-32.27, 37.28, 'LS23')
b.show()
b.save_map('mymap.png', 'My Map')