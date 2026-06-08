#!env python3
"""Read webb format data

A WEBB file contains raw WEBB OBS and pressure gauge data.
This format was used on "Spahr Webb" OBSs from about 1988 to 2000.
Each file contains multiple records and all channels for one instrument,
as well as time stamps, but no instrument response information

Translated from matlab readWebb (1999-2009)
"""
import numpy as np
from copy import copy

from obspy.core import Stream, Trace, Stats
from obspy.core import UTCDateTime

DATA_FORMAT = 'int32'  # necessary to save as STEIM1/2

class WebbDataFile:
    """
    Webb format OBS data (~1991 to 2002)

    Attributes:
        n_channels (int):
        bytes_per_record (int)
        samples_per_record (int):
        samples_per_data_block (list of int):
        sample_rate (float):
        bytes_per_time_block (int):
        time_byte_offset (int):
        bytes_per_sample (int):
        ch_subchannels (list of int):
        ch_names (list of str): each str is "<CHAN>:<LOC>"
        dtype (str): 'u2'
        endian_str (str): '<' for little-endian, '>' for big-endian
        clock_rate (int):
        clock_timezero (UTCDateTime):
        verbose (bool):
    """
    def __init__(self, n_channels, bytes_per_record, samples_per_record,
                 samples_per_data_block, clock_rate, sample_rate,
                 endian, ch_names, bytes_per_time_block=16,
                 time_byte_offset=None, ch_subchannels=None, **kwargs):
        """
        Set up Webb file parameters

        Args:
            n_channels (int): the number of channels
            bytes_per_record (int): the number of bytes per record
            samples_per_record (int):	the number of samples per channel per
                    record
            samples_per_data_block (list of int): Number of samples before
                    and between each header block
            sample_rate (float): the data sampling rate (sps)
            endian (str): 'ieee-be' or '>' (Motorola:STS-DOG),
                          'ieee-le' of '>' (most Webb instruments)
            ch_names (list): the CHAN:LOC code for each channel
            clock_rate (float): the clock rate (for the counter).
                    Usually 256 Hz, but different for some experiments
            bytes_per_time_block (int): Number of bytes in each header block
                    (default=16)
            time_byte_offset (int): Start address (in bytes) of the
                   time byte within time block
            ch_subchannels (list): Number of multiplexed channels
                    within each named channel (default: 1 per channel)
        """
        # Validate parameters
        assert len(ch_names) == n_channels, f'{len(ch_names)=}, {n_channels=}'
        for ch_name in ch_names:
            if isinstance(ch_name, str):
                assert ':' in ch_name, f'{ch_name=}'
            else:
                for name in ch_name:
                    assert ':' in name, f'{name=}'
        assert endian in ('ieee-be', 'ieee-le', '>', '<')
        assert bytes_per_time_block in (4, 8, 16)
        if not len(kwargs) == 0:
            print(f'Extra parameters {list(kwargs.keys())} ignored')

        # Create attributes
        self.bytes_per_record = bytes_per_record
        self.samples_per_record = samples_per_record
        self.samples_per_data_block = samples_per_data_block
        self.clock_rate = clock_rate
        self.sample_rate = sample_rate
        self.n_channels = n_channels
        self.bytes_per_time_block = bytes_per_time_block
        self.time_byte_offset = time_byte_offset
        self.ch_names = ch_names
        self.bytes_per_sample = 2
        self.dtype = 'u2'
        if endian == 'ieee-le' or endian == '<':
            self.endian = '<'
        elif endian == 'ieee-be' or endian == '>':
            self.endian = '>'
        if self.time_byte_offset is None:
            if self.bytes_per_time_block == 4:
                self.time_byte_offset = 0
            elif self.bytes_per_time_block == 8:
                self.time_byte_offset = 2
            elif self.bytes_per_time_block == 16:
                self.time_byte_offset = 4
                
    def __str__(self):
        s = "WebbFile:\n"
        s += f"\tbytes_per_record: {self.bytes_per_record}\n"
        s += f"\tsamples_per_record: {self.samples_per_record}\n"
        s += f"\tsamples_per_data_block: {self.samples_per_data_block}\n"
        s += f"\tclock_rate: {self.clock_rate}\n"
        s += f"\tsample_rate: {self.sample_rate}\n"
        s += f"\tn_channels: {self.n_channels}\n"
        s += f"\tbytes_per_time_block: {self.bytes_per_time_block}\n"
        s += f"\ttime_byte_offset: {self.time_byte_offset}\n"
        s += f"\tbytes_per_sample: {self.bytes_per_sample}\n"
        s += f"\tch_names: {self.ch_names}\n"
        s += f"\tdtype: '{self.dtype}'\n"
        s += f"\tendian: '{self.endian}'"
        return s

    def read(self, fname, station_code, clock_timezero,
             network_code='XX', verbose=False):
        """
        Reads Webb-format OBS data (~1991-2002) into an obspy stream

        For now, reads the whole thing

        Args:
            fname  (str): File name
            station_code (str): Station name
            network_code (str): Network code (default='XX')
            clock_timezero (str): instrument clock time zero in a format that
                obspy.core.UTCDateTime can read
            verbose (bool): Be chatty
            # start (float or UTCDateTime): if a float, start that many seconds
            #    into the file. If an UTCDateTime, start at given time
            # readlen (float): number of seconds to read (0: read to end of file)
        """
        self.station_code = station_code
        self.network_code = network_code
        self.fname = fname
        self.clock_timezero = UTCDateTime(clock_timezero)
        self.verbose = verbose

        print(f'Reading Webb file  {self.fname}, from start to end')
        firstTime = self._read_record_starttime(0)  # Clock drifts not handled
        if verbose is True:
            print(f'firstTime={firstTime.isoformat()}')
        n_recs, n_samps = self._get_file_info()

        readlen = n_samps / self.sample_rate
        print(f'Reading {n_recs} records ({_seconds_string(readlen)})')

        lastAddress = n_samps - 1
        last_rec = n_recs - 1
        stream = self._read_record(0)
        x = range(1, last_rec)
        for recno in x:
            stream += self._read_record(recno)
        return stream

    def _read_record_starttime(self, recNo):
        """
        Return the record start Time

        Args:
            recno (int): record number (first record is 0)
        """
        if self.verbose:
            print('read_webb:in readRecStartTime()')

        with open(self.fname, 'rb') as fid:
            if recNo > 0:
                fid.seek(recNo * self.bytes_per_record, 0)

            # SKIP TO FIRST TIME BLOCK
            skip_samples = self.samples_per_data_block[0]  # skip 1st data block
            skip_bytes = skip_samples * self.n_channels * self.bytes_per_sample
            fid.seek(skip_bytes , 1)
            # READ TIME IN FIRST TIME BLOCK
            # ASSUME Time counter at center of next sample (true in C code)
            timeaddr = skip_samples
            if self.verbose:
                print(f'time address = {timeaddr:8.1f}')
            clock_seconds = self._read_time_block(fid)
            # shift clock ticks back to center of sample #1
            start_clock_secs = clock_seconds - timeaddr / self.sample_rate

            if self.verbose:
                print(f'record {recNo} start clock seconds: {start_clock_secs:.3f}')

        return self.clock_timezero + start_clock_secs

    def _read_time_block(self, fid):
        """
        Reads in time block and returns clock seconds from time zero
        """
        tblc = self.time_byte_offset
        buf = fid.read(self.bytes_per_time_block)
        a = np.frombuffer(buf[tblc:tblc+4], dtype=np.dtype(self.endian + 'u4'))
        clock_counter = a[0]
        clock_secs = clock_counter / self.clock_rate
        if self.verbose:
            print('time: block=0x{}, bytes=0x{:08x}, counts={}, secs={:.3f}'
                  .format(buf.hex(), clock_counter, clock_counter, clock_secs))
        return clock_secs

    def _read_record(self, recno):
        """
        Return one record of data

        Args:
            recno (int): record number (counting from 0)
        """
        with open(self.fname, 'rb') as fid:
            fid.seek(recno*self.bytes_per_record, 0)
            i_spdb = 0  # index of self.samples_per_data_block
            i_data = 0
            firstclock = True
            chdata = np.zeros((self.n_channels, self.samples_per_record), DATA_FORMAT)
            prev_clock_seconds, prev_timeaddr = -1, -1
            while True:
                if (i_data > 0):  # READ A TIME BLOCK
                    # ASSUME time counter corresponds  to the center of the
                    # next sample (consistent with C code)
                    timeaddr = i_data
                    if self.verbose:
                        print('time address = {:8.1f}'.format(timeaddr))
                    clock_seconds = self._read_time_block(fid)
                    if self.verbose:
                        _print_estimated_sample_rate(clock_seconds,
                                                     prev_clock_seconds,
                                                     timeaddr,
                                                     prev_timeaddr)
                        prev_clock_seconds, prev_timeaddr = clock_seconds, timeaddr
                    if firstclock is True:
                        start_clock_secs = clock_seconds - timeaddr/self.sample_rate
                        firstclock = False
                read_samples = self.samples_per_data_block[i_spdb]
                if read_samples + i_data > self.samples_per_record:
                    read_samples = self.samples_per_record - i_data
                read_bytes = read_samples * self.n_channels * self.bytes_per_sample

                buf = fid.read(read_bytes)
                np_data = np.frombuffer(buf, dtype=np.dtype(self.endian + self.dtype))
                a = np.reshape(np_data, (read_samples, self.n_channels)).T
                chdata[:, i_data: i_data + read_samples] = a
                i_data += read_samples
                if i_spdb < len(self.samples_per_data_block) - 1:
                    i_spdb += 1    # Advance to next data block length value
                if i_data >= self.samples_per_record - 1:
                    break
            if self.verbose:
                print(f'{start_clock_secs=:.3f}')

        traces = []
        stats = Stats({'network': self.network_code,
                       'sampling_rate': self.sample_rate,
                       'station': self.station_code,
                       'starttime': self.clock_timezero + start_clock_secs})
        for i, ch_name in zip(range(self.n_channels), self.ch_names):
            if isinstance(ch_name, str):
                stats.channel, stats.location = ch_name.split(':')
                t = Trace(data=chdata[i, :].astype(DATA_FORMAT), header=stats)
                t.data = chdata[i, :].astype(DATA_FORMAT)
                traces.append(t)
            else:
                n_subchannels = len(ch_name)
                ch_len = len(chdata[i,:])
                res = ch_len % n_subchannels
                if res == 0:
                    subchdata = np.reshape(chdata[i, :], (-1, n_subchannels)).T
                else:
                    subchdata = np.reshape(chdata[i, :-res], (-1, n_subchannels)).T
                subch_stats = copy(stats)
                subch_stats.sampling_rate=stats.sampling_rate/n_subchannels
                # Extract subchannel_names and clean up main ch_names list
                for (j, subch_name) in zip(range(n_subchannels), ch_name):
                    subch_stats.channel, subch_stats.location = subch_name.split(':')
                    t = Trace(data=subchdata[j, :].astype(DATA_FORMAT),
                              header=subch_stats)
                    t.data = subchdata[j, :].astype(DATA_FORMAT)
                    traces.append(t)
        stream = Stream(traces=traces)
        return stream

    def _get_file_info(self):
        """Return file length in records and samples"""
        with open(self.fname, 'rb') as fid:
            fid.seek(0, 2)  # Go to the end of the file
            filebytes = fid.tell()
        numRecs = int(np.floor((filebytes-1)/self.bytes_per_record)+1)
        numSamps = numRecs * self.samples_per_record
        return numRecs, numSamps

def _seconds_string(seconds):
    if seconds > 86400:
        days = int(np.floor(seconds/86400))
        s = f"{days}d,"
        seconds -= days*86400
    else:
        s = ''
    hours = int(np.floor(seconds/3600))
    seconds -= hours*3600
    minutes = int(np.floor(seconds/60))
    seconds -= minutes*60
    return s + f"{hours}h,{minutes}m,{seconds}s"


def _print_estimated_sample_rate(clock_seconds, prev_clock_seconds,
                                 time_addr, previous_time_addr):
    if prev_clock_seconds > 0:
        delta_t = clock_seconds-prev_clock_seconds
        delta_i = time_addr-previous_time_addr
        sps = delta_i/delta_t
        print(f'delta = ({delta_i} samples, {delta_t:g} secs) => ({sps:g} samples/sec)')
