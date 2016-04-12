import numpy as np

from scipy.interpolate import UnivariateSpline


class FMdata:
    def __init__(self, lcs, times, n_min_data=100):
        """
        Filling missing values for each light curve using quadratic interpolation
        without smoothing. Note that extrapolation must be strictly prohibited,
        and thus we sync each light curve based on the latest start epoch
        and the earliest end epoch among the all light curves.
        :param lcs: A list of light curves.
        :param times: A list of times.
        :param n_min_data: The minimum number of data points in each
        light curve. If fewer than this, the light curve will be discarded.
        """
        # Type check.
        if type(lcs) != np.ndarray:
            lcs = np.array(lcs)
        if type(times) != np.ndarray:
            times = np.array(times)

        # Dimension check.
        if lcs.shape[0] != times.shape[0]:
            raise RuntimeError('The number of light curves and ' +
                               'the number of times do not match.')

        # Discard light curves having fewer data points than "n_min_data".
        keep_index = []
        for i in range(len(lcs)):
            if len(lcs[i]) >= n_min_data:
                keep_index.append(i)

        # Initialize.
        self.lcs = lcs[keep_index]
        self.times = times[keep_index]

    def run(self):
        # Sync times.
        self._sync_time()

        # Fill missing values.
        self._fill_missing_values()

        return self.filled_lcs, self.synced_epoch

    def _sync_time(self):
        """
        Walk through times of all light curves and create one-dimensional list
        of times that will be used to fill missing values
        for every light curves. In order to prevent extrapolation,
        we chose the latest start epoch and the earliest end epoch
        as the new epoch range.
        """
        # Get all unique epochs and find
        # the latest start epoch and the earliest end epoch.
        all_epoch = []
        latest_start_epoch = -np.inf
        earliest_end_epoch = np.inf
        for t in self.times:
            all_epoch = np.hstack([all_epoch, t])
            if t[0] > latest_start_epoch:
                latest_start_epoch = t[0]
            if t[-1] < earliest_end_epoch:
                earliest_end_epoch = t[-1]

        all_epoch = np.unique(np.ravel(all_epoch))
        all_epoch.sort()

        # Cut epoch.
        start_index = np.searchsorted(all_epoch, latest_start_epoch)
        end_index = np.searchsorted(all_epoch, earliest_end_epoch) + 1
        epoch = all_epoch[start_index:end_index]

        self.synced_epoch = epoch

    def _fill_missing_values(self):
        """
        Fill missing values.
        :return: filled lcs and synced epoch.
        """
        # Fitting.
        filled_lcs = []
        for i in range(len(self.lcs)):
            # Quadratic spline without smoothing.
            spl = UnivariateSpline(self.times[i], self.lcs[i], k=1., s=0.)
            filled_lc = spl(self.synced_epoch)
            filled_lcs.append(filled_lc)

        self.filled_lcs = np.array(filled_lcs)


if __name__ == '__main__':
    # Number of light curves.
    n_lc = 100
    # Number of data points.
    n_data = 200
    # The maximum ratio of missing data points.
    missing_ratio = 0.3

    times = []
    lcs = []
    for i in range(n_lc):
        # The number of missing data points.
        n_missing = int(np.random.rand() * missing_ratio * n_data)
        index = np.arange(0, n_data)
        np.random.shuffle(index)

        t = np.arange(0, n_data)[index[n_missing:]]
        t.sort()
        lc = np.random.randn(n_data)[index[n_missing:]] * 10. + 100.

        times.append(t)
        lcs.append(lc)

    fmdata = FMdata(lcs, times)
    lcs, epoch = fmdata.run()

    from detrend import PDTrend
    pdt = PDTrend(lcs, n_min_member=10, dist_cut=0.5)
    pdt.run()

    import pylab as pl
    pl.plot(epoch, pdt.master_trends[0])
    pl.show()

