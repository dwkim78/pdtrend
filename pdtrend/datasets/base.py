__author__ = 'kim'

"""
Base IO code for all datasets
"""

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
    lcs = pickle.load(bz2.BZ2File(file_path, 'r'))

    return lcs
