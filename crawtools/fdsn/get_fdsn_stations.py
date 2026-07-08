#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download and save stations from IRIS PiLAB

:author: Wayne Crawford (crawford@ipgp.fr), 2019
:requires: obspy>=1.1
"""
import argparse
import sys
import os.path

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

# Global variables
client_address = "IRIS"
formats_suffixes = {'STATIONXML': '.xml',
                    'CSS': '.css',
                    'KML': '.kml',
                    'SACPZ': '.pz',
                    'SHAPEFILE': '.shape',
                    'STATIONTXT': '.txt'}


def main():
    """Main loop"""
    args, kwargs = _get_args()
    client = Client(args.client)
    print('Getting inventory from FDSN server...')
    print(kwargs)
    inv = client.get_stations(**kwargs)
    if args.verbose:
        _print_station_details(inv)
    outfile = _make_outfilename(args)
    if args.plot:
        _plot_stations(inv, args, outfile)
    print(f'Saving to {args.format} file {outfile}')
    inv.write(outfile, format=args.format)


def _get_args():
    """
    Get command line arguments
    """
    formats = [x for x in formats_suffixes.keys()]
    levels = ['network', 'station', 'channel', 'response']
    
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("start_time", type=_to_UTCDateTime,
                        help='Start time (YYYY-MM-DD)')
    parser.add_argument("end_time", type=_to_UTCDateTime,
                        help='End time (YYYY-MM-DD)')
    parser.add_argument("-o", "--outfile", help="Output file name")
    parser.add_argument("-c", "--client", default="IRIS",
                        help='Client (%defaults)s')
    parser.add_argument("-f", "--format", default='STATIONXML', metavar='',
                        choices=formats,
                        help='Output file format.  Allowed values are '
                             +', '.join(formats)
                             + ". (default: %(default)s)")
    parser.add_argument("-b", "--bounds", nargs=4,
                        metavar=('minlat', 'maxlat', 'minlon', 'maxlon'),
                        help="Set latitude/longitude bounds")
    parser.add_argument("-n", "--network", type=str, default=None,
                        help="Limit to the named network")
    parser.add_argument("-d", "--level", default='response', metavar="",
                        choices=levels,
                        help="Limit to given level.  Allowed values are "
                             +", ".join(levels)
                             +". (default: %(default)s)")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Plot stations to file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")
    args = parser.parse_args()

    kwargs = {"starttime": args.start_time,
              "endtime": args.end_time,
              "level": args.level}
    if args.bounds:
        assert args.bounds[0] < args.bounds[1]
        assert args.bounds[2] < args.bounds[3]
        assert args.bounds[0] >= -90.
        assert args.bounds[1] <= 90.
        assert args.bounds[2] >= -180.
        assert args.bounds[3] <= 180.
        kwargs['minlatitude'] = args.bounds[0]
        kwargs['maxlatitude'] = args.bounds[1]
        kwargs['minlongitude'] = args.bounds[2]
        kwargs['maxlongitude'] = args.bounds[3]
    if args.network:
        kwargs['network'] = args.network
    return args, kwargs


def _plot_stations(inv, args, outfile):
    resolution = 'l'
    if args.bounds:
        resolution = 'h'
    plotfile = os.path.splitext(outfile)[0] + '.png'
    print(f'Plotting station map to {plotfile}')
    inv.plot(projection='local', resolution=resolution, size=100,
             outfile=plotfile,  dpi=600, color_per_network=True)


def _make_outfilename(args):
    """Make an output filename"""
    if args.outfile:
        return args.outfile
    fname = 'stations_fdsn'
    if args.bounds:
        fname += '_local'
    fname += '_{}_{}'.format(args.start_time.strftime('%Y%m%d'),
                             args.end_time.strftime('%Y%m%d'))
    if not args.level == 'response':
        fname += f'_{args.level}'
    fname += formats_suffixes.get(args.format, '')
    return fname


def _print_station_details(inv):
    print('Networks = {}'.format(list(set([n.code for n in inv]))))
    print('  Net.Station| Longitude | Latitude   | Start Time       | End_Time')
    print('=============+===========+============+' + 16*'=' + '+' + 16*'=')
    for net in sorted(inv, key=lambda i: i.code):
        for sta in sorted(net, key=lambda i: (i.start_date, i.code)):
            ststr = 'None'
            endstr = 'None'
            if sta.start_date:
                ststr = sta.start_date.strftime('%Y-%m-%dT%H:%M')
            if sta.end_date:
                endstr = sta.end_date.strftime('%Y-%m-%dT%H:%M')
            fmt = ' {:>5s}.{:<5s} | {:<9g} | {:<10g} | {:16s} | {:16s}'
            print(fmt.format(net.code, sta.code, sta.longitude, sta.latitude,
                  ststr, endstr))


def _to_UTCDateTime(str):
    """
    Convert parser input string to UTCDateTime or vomit
    """
    try:
        val = UTCDateTime(str)
    except ValueError:
        print(f'{str} is invalid')
        sys.exit(1)
    else:
        return val


if __name__ == "__main__":
    main()
