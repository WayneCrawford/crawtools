import plot_bathymap as pb
m = pb.plot_basemap(
    0.05, 0.05,
    [-32.4, -32.2, 37.2, 37.4],
    'LuckySISMOMAR_40mtr.grd')
m = pb.plot_station(m, -32.27, 37.25, 'LSVVVVV')
m = pb.plot_station(m, -32.27, 37.28, 'LS23')
pb.save_map('mymap.png', 'My Map', show=True)