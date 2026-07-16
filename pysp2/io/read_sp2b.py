import xarray as xr
import struct
import numpy as np
import platform
import zipfile

from datetime import datetime, timezone


def read_sp2(file_name, debug=False, arm_convention=True):
    """
    Loads a binary SP2 raw data file and returns all of the wave forms
    into an xarray Dataset.

    Parameters
    ----------
    file_name: str
        The name of the .sp2b file to read.
    debug: bool
        Set to true for verbose output.
    arm_convention: bool
        If True, then the file name will follow ARM standard naming conventions.
        If False, then the file name follows the SP2 default naming convention.

    Returns
    -------
    dataset: xarray.Dataset
        The xarray Dataset containing the raw SP2 waveforms for each particle.
    """

    my_data = open(file_name, "rb").read()
    # Get file date from name
    if platform.system() == "Windows":
        split_file_name = file_name.split("\\")
    else:
        split_file_name = file_name.split("/")
    if arm_convention:
        next_split = split_file_name[-1].split(".")
        dt = datetime.strptime(next_split[2], "%Y%m%d")
    else:
        dt = datetime.strptime(split_file_name[-1][0:8], "%Y%m%d")

    if len(my_data) > 0:
        bytepos = 0
        numCols = struct.unpack(">I", my_data[bytepos:bytepos + 4])[0]
        bytepos += 4
        numChannels = struct.unpack(">I", my_data[bytepos:bytepos + 4])[0]
        if debug:
            print(("Loaded file with numCols = {}, numChannels = {}"
                   .format(numCols, numChannels)))

        data_points_per_record = numChannels * numCols

        bytes_per_record = 2 * data_points_per_record
        bytes_not_data_array = 12 + 2 + 28 + 16
        bytes_per_record += bytes_not_data_array
        last_pos = int(bytes_per_record - 1)
        num_spare_cols = struct.unpack(">I", my_data[last_pos - 4:last_pos])[0]
        if debug:
            print("Number of spare columns = %d" % num_spare_cols)

        if num_spare_cols != 0:
            bytes_per_record += num_spare_cols

        numRecords = int(len(my_data) / bytes_per_record)
        totalRows = numChannels * numRecords
        DataWave = np.zeros((totalRows, numCols), dtype='int16')
        Flag = np.zeros(int(totalRows / numChannels), dtype='int16')
        TimeWave = np.zeros(numRecords, dtype='float64')
        Res1 = np.zeros(numRecords, dtype='float32')
        EventIndex = np.zeros(numRecords, dtype='float32')
        TimeDiv10000 = np.zeros(numRecords, dtype='float64')
        TimeRemainder = np.zeros(numRecords, dtype='float64')
        Res5 = np.zeros(numRecords, dtype='float32')
        Res6 = np.zeros(numRecords, dtype='float32')
        Res7 = np.zeros(numRecords, dtype='float64')
        Res8 = np.zeros(numRecords, dtype='float64')
        if num_spare_cols != 0:
            SpareDataArray = np.zeros(numRecords, num_spare_cols)

        arrayFmt = ">"
        for i in range(data_points_per_record):
            arrayFmt += "h"

        for record in range(numRecords):
            dataStartPoint = record * bytes_per_record + 8
            startRow = record * numChannels
            endRow = startRow + numChannels - 1
            the_row = np.array(struct.unpack(
                arrayFmt, my_data[dataStartPoint:dataStartPoint + int(data_points_per_record * 2)]))

            DataWave[startRow:endRow + 1, 0:numCols] = the_row.reshape(
                numCols, numChannels).T
            dataStartPoint += data_points_per_record * 2
            Flag[record] = struct.unpack(">h", my_data[dataStartPoint:dataStartPoint + 2])[0]
            next_floats = struct.unpack(">ffffffff", my_data[dataStartPoint + 2:dataStartPoint + 34])
            TimeWave[record] = next_floats[0]
            Res1[record] = next_floats[1]
            EventIndex[record] = next_floats[2]
            TimeDiv10000[record] = next_floats[3]
            TimeRemainder[record] = next_floats[4]
            Res5[record] = next_floats[5]
            Res6[record] = next_floats[6]
            next_doubles = struct.unpack(">dd", my_data[dataStartPoint + 34:dataStartPoint + 50])
            Res7[record] = next_doubles[0]
            Res8[record] = next_doubles[1]
            dataStartPoint += 50

            if num_spare_cols != 0:
                startRow = (2 * num_spare_cols) * record
                dataStartPoint += bytes_not_data_array - 4
                spareFmt = ">"
                for i in range(num_spare_cols):
                    spareFmt += "f"

                SpareDataArray[record] = np.array(
                    struct.unpack(spareFmt, my_data[dataStartPoint:dataStartPoint+4*num_spare_cols]))

        UTCtime = TimeDiv10000 * 10000 + TimeRemainder
        diff_epoch_1904 = (
            datetime(1970, 1, 1) - datetime(1904, 1, 1)).total_seconds()
        UTCdatetime = np.array([
            datetime.fromtimestamp(
                x - diff_epoch_1904, tz=timezone.utc).replace(tzinfo=None) for x in UTCtime])

        DateTimeWave = (dt - datetime(1904, 1, 1)).total_seconds() + TimeWave

        # Make an xarray dataset for SP2
        Flag = xr.DataArray(Flag, dims={'event_index': EventIndex})
        Res1 = xr.DataArray(Res1, dims={'event_index': EventIndex})
        Res5 = xr.DataArray(Res5, dims={'event_index': EventIndex})
        Res6 = xr.DataArray(Res6, dims={'event_index': EventIndex})
        Res7 = xr.DataArray(Res7, dims={'event_index': EventIndex})
        Res8 = xr.DataArray(Res8, dims={'event_index': EventIndex})
        Time = xr.DataArray(UTCdatetime, dims={'event_index': EventIndex})
        EventInd = xr.DataArray(EventIndex, dims={'event_index': EventIndex})
        DateTimeWaveUTC = xr.DataArray(UTCtime, dims={'event_index': EventIndex})
        DateTimeWave = xr.DataArray(DateTimeWave, dims={'event_index': EventIndex})
        TimeWave = xr.DataArray(TimeWave, dims={'event_index': EventIndex})
        my_ds = xr.Dataset({'time': Time, 'Flag': Flag, 'Res1': Res1, 'Res5': Res5,
                            'Res6': Res6, 'Res7': Res7, 'Res8': Res8, 'EventIndex': EventInd,
                            'DateTimeWaveUTC': DateTimeWaveUTC, 'TimeWave': TimeWave,
                            'DateTimeWave': DateTimeWave})

        for i in range(numChannels):
            temp_array = np.zeros((numRecords, numCols), dtype='int')
            for j in range(numRecords):
                k = i + j*numChannels
                temp_array[j] = DataWave[k]
            my_ds['Data_ch' + str(i)] = xr.DataArray(
                temp_array, dims={'event_index': EventIndex, 'columns': np.arange(0, 100, 1)})
        del my_data
        del DataWave
        return my_ds
    else:
        return None


def read_sp2xr(file_name, debug=False):
    """
    Loads a binary SP2-XR raw data file and returns all of the waveforms
    into an xarray Dataset.

    The SP2-XR file format differs from the original SP2:
    each record carries its own (rows, cols) header instead of a single
    file-level header, the waveform samples are stored as 4-byte integers
    (rather than 2-byte) and only two detector channels are written
    (ch0 = scattering, ch1 = incandescence). The time encoding is the
    same as the SP2 (TimeDiv10000 * 10000 + TimeRemainder, in seconds
    since 1904-01-01).

    Parameters
    ----------
    file_name: str
        The name of the .sp2b file (or a .zip containing the .sp2b) to read.
    debug: bool
        Set to true for verbose output.

    Returns
    -------
    dataset: xarray.Dataset
        The xarray Dataset containing the raw SP2-XR waveforms for each
        particle. Variable names mirror read_sp2() (Data_ch0, Data_ch1,
        Flag, EventIndex, time, ...).
    """

    if file_name.lower().endswith(".zip"):
        with zipfile.ZipFile(file_name, "r") as zf:
            inner = zf.namelist()[0]
            my_data = zf.read(inner)
    else:
        my_data = open(file_name, "rb").read()

    if len(my_data) == 0:
        return None

    # Pre-scan: locate the start byte of every record and verify shape consistency
    record_starts = []
    pos = 0
    numSamples = None
    numChannels = None

    while pos + 8 <= len(my_data):
        rows, cols = struct.unpack(">2I", my_data[pos:pos + 8])
        # rows == 0 marks an explicit end-of-data marker
        if rows == 0:
            break

        if numSamples is None:
            numSamples, numChannels = rows, cols
        elif (rows, cols) != (numSamples, numChannels):
            raise ValueError(
                "SP2-XR records with varying waveform dimensions are not "
                "supported (record %d: (%d, %d) vs first (%d, %d))." %
                (len(record_starts), rows, cols, numSamples, numChannels))

        wave_bytes = rows * cols * 4
        # waveform + Flag(2) + 7 floats(28) + 2 doubles(16) + array_2_dim(4)
        meta_end = pos + 8 + wave_bytes + 2 + 28 + 16 + 4
        if meta_end > len(my_data):
            break
        a2_dim = struct.unpack(">I", my_data[meta_end - 4:meta_end])[0]
        rec_end = meta_end + a2_dim * 4
        if rec_end > len(my_data):
            break

        record_starts.append(pos)
        pos = rec_end

    numRecords = len(record_starts)
    if numRecords == 0:
        return None

    if debug:
        print("Loaded file with numRecords = %d, numChannels = %d, "
              "numSamples = %d" % (numRecords, numChannels, numSamples))

    DataWave = np.zeros((numRecords, numChannels, numSamples), dtype='int32')
    Flag = np.zeros(numRecords, dtype='int32')
    TimeWave = np.zeros(numRecords, dtype='float32')
    Res1 = np.zeros(numRecords, dtype='float32')
    EventIndex = np.zeros(numRecords, dtype='float32')
    TimeDiv10000 = np.zeros(numRecords, dtype='float32')
    TimeRemainder = np.zeros(numRecords, dtype='float32')
    Res5 = np.zeros(numRecords, dtype='float32')
    Res6 = np.zeros(numRecords, dtype='float32')
    Res7 = np.zeros(numRecords, dtype='float64')
    Res8 = np.zeros(numRecords, dtype='float64')

    wave_bytes = numSamples * numChannels * 4
    wave_fmt = ">" + "i" * (numSamples * numChannels)

    for i, start in enumerate(record_starts):
        # Skip the per-record (rows, cols) header
        p = start + 8

        # Waveform: interleaved (ch0[0], ch1[0], ch0[1], ch1[1], ...)
        raw = np.array(struct.unpack(wave_fmt, my_data[p:p + wave_bytes]),
                       dtype='int32')
        # Reshape (numSamples, numChannels) then transpose to (numChannels, numSamples)
        DataWave[i] = raw.reshape(numSamples, numChannels).T
        p += wave_bytes

        Flag[i] = struct.unpack(">H", my_data[p:p + 2])[0]
        p += 2

        floats = struct.unpack(">7f", my_data[p:p + 28])
        TimeWave[i] = floats[0]
        Res1[i] = floats[1]
        EventIndex[i] = floats[2]
        TimeDiv10000[i] = floats[3]
        TimeRemainder[i] = floats[4]
        Res5[i] = floats[5]
        Res6[i] = floats[6]
        p += 28

        doubles = struct.unpack(">2d", my_data[p:p + 16])
        Res7[i] = doubles[0]
        Res8[i] = doubles[1]

    # Reconstruct UTC datetime exactly like read_sp2()
    UTCtime = TimeDiv10000.astype('float64') * 10000 + TimeRemainder.astype('float64')
    diff_epoch_1904 = (
        datetime(1970, 1, 1) - datetime(1904, 1, 1)).total_seconds()
    UTCdatetime = np.array([
        datetime.fromtimestamp(
            x - diff_epoch_1904, tz=timezone.utc).replace(tzinfo=None)
        for x in UTCtime])

    # Build xarray Dataset
    Flag_da = xr.DataArray(Flag, dims={'event_index': EventIndex})
    Res1_da = xr.DataArray(Res1, dims={'event_index': EventIndex})
    Res5_da = xr.DataArray(Res5, dims={'event_index': EventIndex})
    Res6_da = xr.DataArray(Res6, dims={'event_index': EventIndex})
    Res7_da = xr.DataArray(Res7, dims={'event_index': EventIndex})
    Res8_da = xr.DataArray(Res8, dims={'event_index': EventIndex})
    Time = xr.DataArray(UTCdatetime, dims={'event_index': EventIndex})
    EventInd = xr.DataArray(EventIndex, dims={'event_index': EventIndex})
    DateTimeWaveUTC = xr.DataArray(UTCtime, dims={'event_index': EventIndex})
    TimeWave_da = xr.DataArray(TimeWave, dims={'event_index': EventIndex})

    my_ds = xr.Dataset({
        'time': Time, 'Flag': Flag_da,
        'Res1': Res1_da, 'Res5': Res5_da, 'Res6': Res6_da,
        'Res7': Res7_da, 'Res8': Res8_da,
        'EventIndex': EventInd,
        'DateTimeWaveUTC': DateTimeWaveUTC,
        'TimeWave': TimeWave_da})

    columns = np.arange(numSamples)
    for ch in range(numChannels):
        my_ds['Data_ch' + str(ch)] = xr.DataArray(
            DataWave[:, ch, :],
            dims={'event_index': EventIndex, 'columns': columns})

    return my_ds
