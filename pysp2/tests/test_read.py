import pysp2


def test_read_sp2b():
    my_sp2 = pysp2.io.read_sp2(pysp2.testing.EXAMPLE_SP2B)
    assert my_sp2.dims['event_index'] == 5877


def test_read_config():
    my_config = pysp2.io.read_config(pysp2.testing.EXAMPLE_INI)
    assert int(my_config['Acquisition']['1 of Every']) == 1


def test_read_hk():
    my_hk = pysp2.io.read_hk_file(pysp2.testing.EXAMPLE_HK)
    assert my_hk['Duty Cycle'].max() == 720.
    assert my_hk['Duty Cycle'].min() == 0.
