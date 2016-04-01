import numpy as np
import pylab as pl

from pdtrend.detrend import PDTrend
from pdtrend.datasets.base import load_lightcurve_set
from pdtrend.utils.logger import Logger


def test():
    logger = Logger().getLogger()

    logger.info('Loading the light curve set.')
    lcs = load_lightcurve_set()
    logger.info('\tThe number of light curves is %d.' % len(lcs))

    logger.info('Initializing pdtrend.')
    pdt = PDTrend(lcs)

    # We can run all the following routines using "pdtrend.run()",
    # but, here we do it individually to print log messages.
    logger.info('Calculating the distance matrix.')
    pdt._calculate_distance_matrix()
    logger.info('Searching for clusters using Birch.')
    pdt._find_clusters()
    logger.info('Filtering the clusters.')
    pdt._filter_clusters()
    logger.info('Building master trends.')
    pdt._build_master_trends()

    logger.info('Detrending one light curves using the master trends.')
    detrended = pdt.detrend(lcs[1])

    logger.info('Ploting results.')
    pl.figure(figsize=(14, 8))
    pl.subplot(211)
    pl.plot(lcs[1], 'b.', label='Raw')
    pl.text(8955, 70000, r'$\sigma$: %.1f' % np.std(lcs[1]), fontsize=15,
            va='center', bbox=dict(boxstyle='round', ec='w', fc='g', alpha=0.3))
    pl.ylabel('Flux')
    pl.grid()
    pl.legend(numpoints=1)

    pl.subplot(212)
    pl.plot(detrended, 'b.', label='Detrended')
    pl.text(8955, 70000, r'$\sigma$: %.1f' % np.std(detrended), fontsize=15,
            va='center', bbox=dict(boxstyle='round', ec='w', fc='g', alpha=0.3))
    pl.ylabel('Flux')
    pl.xlabel('Time index')
    pl.grid()
    pl.legend(numpoints=1)
    pl.show()

    logger.info('Done.')
    logger.handlers = []


if __name__ == '__main__':
    test()