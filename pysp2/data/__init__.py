import os
import pandas as pd

SCAT_TABLE_PATH = os.path.join(os.path.dirname(__file__), 'sp2MieScatteringTableoffcoin.v1.xls')
CalibrationTable = pd.read_excel(SCAT_TABLE_PATH, header=1)