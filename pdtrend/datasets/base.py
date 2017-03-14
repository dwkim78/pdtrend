"""
Base IO code for all datasets
"""

import sys

from os.path import dirname
from os.path import join


def load_lightcurve_set():
    """
    Return the set of light curves for testing pdtrend.

    Returns
    -------
    lcs : numpy.ndarray
        An array of light curves.
    """

    import bz2

    try:
        import cPickle as pickle
    except:
        import pickle

    module_path = dirname(__file__)

    # The light curves are bzipped and pickled.
    file_path = join(module_path, 'lightcurves/lc.pbz2')
    # For Python 3.
    if sys.version_info.major >= 3:
        lcs = pickle.load(bz2.BZ2File(file_path, 'r'), encoding='bytes')
    # For Python 2.
    else:
        lcs = pickle.load(bz2.BZ2File(file_path, 'r'))

    return lcs
