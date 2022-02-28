============================================
Information contained in particle wave files
============================================

.. list-table:: Variables extracted from .sp2b files
   :width: 50 50
   :header-rows: 1

   * - Variable name
     - Variable description
   * - time
     - timestamp in seconds since the epoch
   * - Flag
     - Internal flag
   * - Res1
     - Res1
   * - Res5
     - Res5
   * - Res6
     - Res6
   * - Res7
     - Res7
   * - Res8
     - Res8
   * - EventIndex
     - The index of the particle event
   * - DateTimeWaveUTC
     - The number of seconds since midnight on Jan 1, 1904.
   * - TimeWave
     - Time stamp in seconds past midnight
   * - DateTimeWave
     - Local time stamp in seconds past midnight on Jan 1, 1904
   * - Data_ch#
     - The waveform in channel #. These numbers range from -32767 to 32767.


.. list-table:: Particle statistics
   :width: 50 50
   :header-rows: 1
   * - Variable name
     - Variable description
   * - DateTimeWaveUTC
     - The number of seconds since midnight on Jan 1, 1904.
   * - TimeWave
     - Time stamp in seconds past midnight
   * - DateTimeWave
     - Local time stamp in seconds past midnight on Jan 1, 1904
   * - Base_ch#
     - The base value of the particle spectrum in channel #.
   * - PkHt_ch#
     - The height of the peak in channel #, after correcting for baseline.
   * - PkPos_ch#
     - The location of the peak (in peak indices, or 0.2 microsecond divisions)
   * - PkStart_ch#
     - The start of the spectrum in peak indices
   * - PkEnd_ch#
     - The end of the spectrum in peak indicies
   * - PkHalfRise_ch#
     - Position of the half peak, on the rising side
   * - PkHalfDecay_ch#
     - Position of the half peak, on the falling side
   * - FtAmp_ch#
     - Gaussian fit amplitude of peak (channels 0, 4 only)
   * - FtPos_ch#
     - Position of Gaussian fit of the peak (channel 0, 4 only)
   * - PkFWHM_ch#
     - The half width of the Gaussian fit of the peak
   * - IncanRatioch5ch6
     - Incandescent temperature ratio, Channel 5/Channel 6
   * - IncanPkOffsetch5ch6
     - Peak Position ch5 - Peak Position ch6
   * - IncanRatioch1ch2
     - Incandescent temperature ratio, Channel 5/Channel 6
   * - IncanPkOffsetch1ch2
     - Peak Position ch5 - Peak Position ch6
