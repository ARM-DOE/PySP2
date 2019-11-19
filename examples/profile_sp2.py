import cProfile
import pysp2

if __name__ == "__main__":
    my_binary = pysp2.io.read_sp2('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110x001.sp2b',
                                  debug=True)
    my_config = pysp2.io.read_config('/nfs/gce/projects/digr/SP2Data4ANL/20181110/20181110114046.ini')
    cProfile.run("pysp2.util.gaussian_fit(my_binary, my_config, parallel=False)")