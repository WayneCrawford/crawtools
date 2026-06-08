# crawtools

Various codes and tools

## Modules:

- `shape_selection`: create map shapes to select/plot events
- `map_cross_section`: create a cross-section map
- `bathy_map`: plot a bathymetric map
    - `BathyMap`: base class
- `instruments`: Simple modeling of instruments (sensor+preamp+logger)


## Other

- `driftcorrect`: codes for clock drift correction:
    - `driftplot`: plot residuals in a QuakeML file as a function of time
    - `driftcorrect`: apply a drift correction one station in a QuakeML file
- `FDSN`: routines to get data/information using FDSN web services.  Only has get_stations for now

Documentation:
https://crawtools.readthedocs.io/