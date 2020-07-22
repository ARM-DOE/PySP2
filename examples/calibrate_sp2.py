import pysp2
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

calibration_path = 'C:\\Users\\rjackson\\Documents\\SP2Data4ANL\\AquaDag_PSL_cal\\pyprocessing\\'
my_cal_data = pysp2.io.read_calibration(calibration_path)
calibration_path_igor = 'C:\\Users\\rjackson\\Documents\\SP2Data4ANL\\AquaDag_PSL_cal\\PAPI_Generated_Numbers\\'
my_cal_data_igor = pysp2.io.read_calibration(calibration_path_igor)
print(my_cal_data.keys())


scat_350 = my_cal_data['scat_350'] # .where(my_cal_data['scat_350'].ScatRejectKey == 0)
scat_350_igor = my_cal_data_igor['scat_350'] # .where(my_cal_data_igor['scat_350'].ScatRejectKey == 0)

print(scat_350)
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
#p = figure(tools="pan, box_zoom, reset, save")
bin_ranges = np.linspace(0, 65535, 10000)

ScatRejectKey = scat_350.ScatRejectKey.values
ScatRejectKey_igor = scat_350_igor.ScatRejectKey.values
ScatRejectKey[ScatRejectKey > 7] = 0
ax.hist(ScatRejectKey, bins=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]) - 0.5, alpha=0.5, label='Python')
ax.hist(ScatRejectKey_igor, bins=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]) - 0.5, alpha=0.5, label='Igor')
ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
ax.legend()
ax.set_xticklabels(['Accepted', 'Invalid', 'Peak\n Low', 'Peak\n High',
                    'Peak\n Narrow', 'Peak\n Wide', 'Peak\n Close', 'Peak\n Far',
                    'Min\n Peak\n Ratio', 'Max\n Peak\n Ratio', 'Peak\n Offset'])
ax.set_xlabel('Rejection reason')
ax.set_ylabel('Count')

plt.show(fig)
