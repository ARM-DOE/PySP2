"""
These are the sample files used for testing PySP2.
"""
import os
from pathlib import Path
from arm_test_data import DATASETS

DATA_PATH = Path(__file__).resolve().parent / "data"


def _sample_file(local_name: str) -> str:
    local_path = DATA_PATH / local_name
    if local_path.exists():
        return str(local_path)
    return DATASETS.fetch(local_name)

#DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')
EXAMPLE_SP2B = _sample_file(
    'mosaossp2M1.00.20191216.130601.raw.20191216x193.sp2b')
EXAMPLE_INI = _sample_file(
    'mosaossp2M1.00.20191216.000601.raw.20191216000000.ini')
EXAMPLE_HK = _sample_file(
    'mosaossp2auxM1.00.20191217.010801.raw.20191216000000.hk')

print("Fetching files from ARM test data repository...")
print("DATASETS:", DATASETS.registry_files)
EXAMPLE_SP2B_PSL = _sample_file("20230721x002.sp2b")
EXAMPLE_INI_PSL = _sample_file("20230721121710.ini")
EXAMPLE_HK_PSL = _sample_file("20230721121711.hk")

SP2XR_DATA_PATH = os.path.join(DATA_PATH, 'SP2XR')
EXAMPLE_SP2XR_SP2B = os.path.join(
    SP2XR_DATA_PATH, 'SP2XR_sp2b_20190508172218_x0001.zip')
EXAMPLE_SP2XR_HK = os.path.join(
    SP2XR_DATA_PATH, 'SP2XR_hk_20190508172218_x0001.zip')
EXAMPLE_SP2XR_INI = os.path.join(
    SP2XR_DATA_PATH, '20190508172218 Calibration 20181005.ini')
EXAMPLE_SP2XR_PBP = os.path.join(
    SP2XR_DATA_PATH, 'SP2XR_PbP_20190508172218_x0001.zip')