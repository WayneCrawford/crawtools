#!/usr/bin/env/python3
"""
Plot EQ catalog time residuals as a function of time
"""
import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path

from obspy import read_events
import matplotlib.pyplot as plt

def _get_stations(cat):
    stations = []
    for event in cat:
        for pick in event.picks:
            if pick.waveform_id.station_code not in stations:
                stations.append(pick.waveform_id.station_code)
    return stations

def _get_pick(pick, origin, phase_regex):
    time, tres = None, None
    if re.match(phase_regex, pick.phase_hint):
        pick_id = pick.resource_id
        for arrival in origin.arrivals:
            if arrival.pick_id == pick_id:
                time = origin.time._get_datetime()
                tres = arrival.time_residual
    return time, tres

def driftplot(catalog_file, stations=None, phases=[r'[pP].*', r'[sS].*']):
    """
    Plots time residuals as a function of time
    
    :param catalog_file: name of QuakeML file with events
    :param stations: list of stations to limit to
    :param phases: regexps of phases to plot (up to 2)
    """
    catalog = read_events(catalog_file)
    if stations is None:
        stations = _get_stations(catalog)
        print(f'found {stations=}')
    
    fig, axs = plt.subplots(len(stations), 1, sharex=True)
    colors='brgmcyk'
    for station, ax in zip(stations, axs):
        times, tress = [], []
        for phase in phases:
            times.append([])
            tress.append([])
        for event in catalog:
            for pick in event.picks:
                if pick.waveform_id.station_code == station:
                    for i in range(len(phases)):
                        time, tres = _get_pick(pick, event.origins[0], phases[i])
                        if time is not None:
                            times[i].append(time)
                            tress[i].append(tres)
        for i in range(len(phases)):
            ax.scatter(times[i], tress[i], s = 2.5, edgecolor = colors[i],
                       label=phases[i])
        ax.axhline(y=0, color='k', linestyle='-')
        ax.set_ylabel(station)
        ax.set_ylim(-1, 1)
    axs[0].set_title("Time residuals")
    axs[0].legend()
    axs[-1].set_xlabel("Origin Time")
    plt.savefig(f"{catalog_file}_drift.png", format='png')
    plt.show()
        
def main():
    parser = ArgumentParser(prog="driftplot", description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("catalog_file", help="Catalog file (QuakeML format)")
    parser.add_argument("--stations", default=None, nargs='+',
                        help="List of stations to limit to")
    parser.add_argument("--p1", default=r'[pP].*',
                        help="First phase to plot (regex)")
    parser.add_argument("--p2", default=r'[sS].*',
                        help="Second phase to plot (regex)")
    args = parser.parse_args()
    
    assert Path(args.catalog_file).is_file()

    driftplot(args.catalog_file, stations=args.stations, phases = [args.p1, args.p2])

if __name__ == "__main__":
    main()

