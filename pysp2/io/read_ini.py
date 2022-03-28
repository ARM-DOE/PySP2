import configparser


def read_config(file):
    """
    This loads the INI file into an easy to use ConfigParser() object.
    This configuration structure can then be made an input to the
    processing functions.

    Parameters
    ----------
    file: str
        The name of the INI file to load

    Returns
    -------
    config: configparser.ConfigParser
        The ConfigParser object containing the structure of the INI file
    """

    my_parser = configparser.ConfigParser()
    my_parser.read(file)

    return my_parser
