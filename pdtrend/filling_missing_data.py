import numpy as np
import warnings

from scipy.interpolate import UnivariateSpline


class FMdata:
    """
    Filling missing values using quadratic interpolation.

    Extrapolation is strictly prohibited, and thus we sync
    each light curve based on the latest start epoch and the earliest end
    epoch among the all light curves.

    Parameters
    ----------
    lcs : array_like
        An array containing a set of light curves. The number of data points
        in each light curve does not need to be identical.
    times : array_like
        An array of observation epochs of the corresponding light curve.
    weights : array_like, optional
        A list of weights. Default is None. Weights are not
        used to fill missing values. FMdata will just discard weights based
        on "n_min_data".
    n_min_data : int, optional
        The minimum number of data points in each light curve.
        If fewer than this, the light curve will be discarded.
    """
    def __init__(self, lcs, times, weights=None, n_min_data=100):
        # Type check.
        if type(lcs) != np.ndarray:
            lcs = np.array(lcs)
        if type(times) != np.ndarray:
            times = np.array(times)
        if weights is not None:
            if type(weights) != np.ndarray:
                weights = np.array(weights)

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
        if weights is not None:
            self.weights = weights[keep_index]
        else:
            self.weights = weights

    def run(self):
        """
        Fill missing values and returns results.

        Returns
        -------
        filled_lcs: numpy.ndarray
            An array of light curves after filling missing values.
        synced_epoch: numpy.ndarray
            Epoch synced between all light curves.
        """
        # Sync times.
        self._sync_time()

        # Fill missing values.
        self._fill_missing_values()
        if self.weights is None:
            return self.filled_lcs, self.synced_epoch
        else:
            return self.filled_lcs, self.synced_epoch, self.weights

    def _sync_time(self):
        """
        Walk through times of all light curves and create one-dimensional
        list of times that will be used to fill missing values for every light
        curves. In order to prevent extrapolation, we chose the latest start
        epoch and the earliest end epoch as the new epoch range.
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
        """Fill missing values for each light curve."""
        filled_lcs = []
        filled_indices = []
        for i in range(len(self.lcs)):
            # Check if the length is same. Print warning but not break the loop.
            # Thus, need to discard the corresponding weight as well.
            if len(self.times[i]) != len(self.lcs[i]):
                warnings.warn(('The number of data points of the %dth ' +
                               'light curve is not matched with the one ' +
                               'of the list of times, and thus is discarded.')
                              % (i + 1))
                continue

            # Quadratic spline without smoothing.
            spl = UnivariateSpline(self.times[i], self.lcs[i], k=1., s=0.)
            filled_lc = spl(self.synced_epoch)
            filled_lcs.append(filled_lc)

            # Indices for the new weight list.
            filled_indices.append(i)

        self.filled_lcs = np.array(filled_lcs)
        if self.weights is not None:
            self.weights = self.weights[filled_indices]
