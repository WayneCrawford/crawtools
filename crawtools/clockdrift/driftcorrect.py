#!/usr/bin/env/python3
"""
Apply a linear drift correction to one station's arrival picks

Correction is zero at starttime, offset at endtime and linear in between
"""
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path

from obspy import read_events, UTCDateTime

def driftcorrect(catalog_file, station, starttime, endtime, offset,
                 out_filename='out_catalog.qml'):
    """
    correct time residuals on one station
    
    :param catalog_file: name of QuakeML file with events
    :param station: station
    :param starttime: sychronization start time (assumed perfectly synchronized)
    :param endtime: synchronization end time
    :param offset: seconds to add to pick times at endtime
    :param out_filename: output catalog name
    """
    catalog = read_events(catalog_file)
    starttime = UTCDateTime(starttime)
    endtime = UTCDateTime(endtime)
    assert endtime > starttime
    sync_span = endtime - starttime
    print(f'{offset=}, {sync_span/86400=}')
    
    firstpick = endtime
    lastpick = starttime
    for event in catalog:
        for pick in event.picks:
            if pick.waveform_id.station_code == station:
                firstpick = min(pick.time, firstpick)
                lastpick = max(pick.time, lastpick)
                pick.time += offset * (pick.time - starttime) / sync_span
    print(f'first, last pick times were {firstpick}, {lastpick}')
    assert firstpick > starttime
    assert lastpick < endtime
    catalog.write(out_filename, format='QUAKEML') 
    print(f'Wrote to {out_filename}')     

def main():
    parser = ArgumentParser(prog="driftcorrect", description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("catalog_file", help="Catalog file (QuakeML format)")
    parser.add_argument("station", help="Station to process")
    parser.add_argument("offset", type=float, help="offset to add to picks at end time")
    parser.add_argument("starttime",
                        help="synchronization start time (when OBS was "
                             "perfectly synchronized)")
    parser.add_argument("endtime", help="synchronization end time")
    parser.add_argument("--outfile", default='out_catalog.qml',
                        help="output catalog filename")

    args = parser.parse_args()
    
    assert Path(args.catalog_file).is_file()

    driftcorrect(args.catalog_file, station=args.station,
                 starttime=UTCDateTime(args.starttime),
                 endtime=UTCDateTime(args.endtime),
                 offset=args.offset, out_filename=args.outfile)

if __name__ == "__main__":
    main()

