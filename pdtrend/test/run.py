import os

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

    np.random.seed(1024)
    # Create random weights.
    #weights = np.random.rand(lcs.shape[0])
    # Same weights.
    weights = np.ones(lcs.shape[0])

    # X and Y coordinates
    xy_coords = np.random.rand(lcs.shape[0], 2) * 1000.

    logger.info('Initializing pdtrend.')
    pdt = PDTrend(lcs, weights=weights, xy_coords=xy_coords,
                  n_min_member=10, dist_cut=0.45)

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

    logger.info('Detrending a light curve using the master trends.')
    detrended = pdt.detrend(lcs[1])

    logger.info('Plotting results.')
    # Creating an output folder.
    output_path = './outputs/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Plotting spatial distribution.
    pdt.plot_spatial('%s/spatial.png' % output_path)

    # Plotting the master trends.
    pl.figure(figsize=(10, 4))
    for i in range(len(pdt.master_trends)):
        pl.subplot(len(pdt.master_trends), 1, i + 1)
        pl.plot(pdt.master_trends[i], 'b.', label='Master trend %d' % (i + 1))
        pl.ylabel('Normalized flux')
        pl.xlabel('Time index')
        pl.legend(numpoints=1, loc='lower right')
        pl.grid()
    pl.tight_layout()
    pl.savefig('%s/master_trends.png' % output_path)

    # Plotting a detrended result of one light curve.
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
    pl.savefig('%s/detrended.png' % output_path)

    logger.info('Done.')
    logger.handlers = []


if __name__ == '__main__':
    test()