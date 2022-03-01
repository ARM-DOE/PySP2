import xarray as xr
import struct
import numpy as np
import platform

from datetime import datetime


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
            datetime.fromtimestamp(x - diff_epoch_1904) for x in UTCtime])

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
