import numpy as np

from sklearn.cluster import Birch
from scipy.optimize import least_squares


class PDTrend:
    """
    Detrending systematic trends in a set of light curves.

    PDTrend determines master trends using Pearson correlation coefficients
    between each light curve and using machine learning, and then provide a
    function to detrend each light curve using the determined master trends.
    For details, see Kim et al. 2009 or visit the GitHub repository
    (https://goo.gl/uRAXfr).

    Parameters
    ----------
    lcs : array_like
        A list of light curves. It must be the shape of N x M,
        where N is the number of light curves and M is the number of
        data points. Note that M must be same for all N light curves.
        pdtrend assumes that these M data points are synced in time.
        Thus pdtrend does not use any time information.
    weights : array_like, optional
        A list of weights for the corresponding light curves.
        It is used only when constructing master trends. It must contain N
        elements. N is the number of the light curves. See the "lcs" parameter.
        Default is None, so the identical weights for all light curves.
    xy_coords : array_like, optional
        X and Y coordinates of stars of the light curves.
        It must contain Nx2 elements, where N is the number of the light curves.
        If the coordinates are given, the function "plot_spatial" can be called
        after constructing master trends, which will plot the spatial
        distribution of the master trends.
    n_min_member : int, optional
        The minimum number of members in each cluster. Default is 10.
    dist_cut : float, optional
        Distance cut to filter clusters found using the Birch clustering
        algorithm. If the median distance
        between members in a cluster is larger than the cut,
        the cluster is discarded. Must be between [0, 1]. Default is 0.45.
    branching_factor : int, optional
        Branching factor for the Birch clustering. Default is 50.
    threshold : float, optional
        Threshold for the Birch clustering. Default is 0.5.
    """
    def __init__(self, lcs, weights=None, xy_coords=None,
                 n_min_member=10, dist_cut=0.45,
                 branching_factor=50, threshold=0.5):
        # Convert the light curve set to numpy array.
        if type(lcs) != np.ndarray:
            lcs = np.array(lcs)

        # Sanity check.
        if len(lcs.shape) != 2:
            raise RuntimeError('lcs must be a 2-dimensional array.')

        if lcs.shape[0] < n_min_member:
            raise RuntimeError('The number of light curves in lcs ' +
                               'is fewer than the n_min_member.')
        if weights is not None:
            if type(weights) != np.ndarray:
                weights = np.array(weights)

            if lcs.shape[0] != weights.shape[0]:
                raise RuntimeError('Shapes of lcs and weights do not match.')
        else:
            # Same weights for the all light curves.
            weights = np.ones(lcs.shape[0])

        if xy_coords is not None:
            if type(xy_coords) != np.ndarray:
                xy_coords = np.array(xy_coords)

            if lcs.shape[0] != xy_coords.shape[0]:
                raise RuntimeError('Shapes of lcs and xy_coords do not match.')

        # Set parameters.
        self.lcs = lcs
        self.weights = weights
        self.xy_coords = xy_coords
        self.n_min_member = n_min_member
        self.dist_cut = dist_cut

        # Initialize.
        self.corr_matrix = None
        self.dist_matrix = None
        self.birch = None
        self.branching_factor = branching_factor
        self.threshold = threshold

    def _calculate_distance_matrix(self):
        """
        Calculate a distance matrix, which is defined as:
        (1. - correlation_matrix) / 2.
        """
        corr_matrix = np.corrcoef(self.lcs)
        dist_matrix = (1. - corr_matrix) / 2.

        self.corr_matrix = corr_matrix
        self.dist_matrix = dist_matrix

    def _find_clusters(self):
        """Find clusters using Birch and the distance matrix."""
        # TODO: Need to test with different threshold.
        # Need to test with multiple dataset having trends.
        # Branching factor is fine.
        birch = Birch(branching_factor=self.branching_factor,
                      threshold=self.threshold,
                      n_clusters=None).fit(self.dist_matrix)

        self.birch = birch

    def _filter_clusters(self):
        """
        Discard a cluster if 1) it has less than "n_min_member" members, or
        2) median distance between each member is larger than "dist_cut".
        """
        unique_labels = set(self.birch.labels_)
        _filtered_labels = []

        for label in unique_labels:
            index = [i for i in range(len(self.birch.labels_)) if
                     self.birch.labels_[i] == label]
            # The number of members in the given cluster.
            if len(index) < self.n_min_member:
                continue

            dist_list = []
            for i in range(len(index) - 1):
                for j in range(i + 1, len(index)):
                    dist_list.append(self.dist_matrix[index[i], index[j]])

            # Median distance check.
            if np.median(dist_list) <= self.dist_cut:
                _filtered_labels.append(label)

        self._filtered_labels = _filtered_labels

        # Check how many clusters are left.
        if len(self._filtered_labels) == 0:
            raise RuntimeWarning(
                'No clusters were found. ' +
                'Adjust input parameters and try again. ' +
                'For instance, decrease "n_min_member" or ' +
                'increase "dist_cut". For details, ' +
                'visit https://github.com/dwkim78/pdtrend'
            )

    def _build_master_trends(self):
        """Build master trends using the filtered clusters."""
        master_trends_indices = []
        master_trends = []
        for label in self._filtered_labels:
            index = [i for i in range(len(self.birch.labels_)) if
                     self.birch.labels_[i] == label]
            master_trends_indices.append(index)

            trends = []
            weights_sum = 0
            for i in range(len(index)):
                # Normalization.
                med_lc = np.median(self.lcs[index[i]])
                normed = (self.lcs[index[i]] - med_lc) / med_lc

                # Weights.
                weights = self.weights[index[i]]
                normed *= weights
                weights_sum += weights

                trends.append(normed)

            # Construct a master trends using the normalized and weighted
            # light curves from the above for loop.
            master_trends.append(np.sum(trends, axis=0) / weights_sum)

        self.master_trends_indices = master_trends_indices
        self.master_trends = master_trends

    def run(self):
        """Run pdtrend pipeline."""
        if self.corr_matrix is None or self.dist_matrix is None:
            self._calculate_distance_matrix()
            self._find_clusters()
        else:
            # For the safety, just in case.
            # If corr_matrix and dist_matrix is not None,
            # birch must be not None either.
            if self.birch is None:
                self._find_clusters()

        self._filter_clusters()
        self._build_master_trends()

    def _func_trends(self, p, x):
        """Return sum of the trends."""
        return np.sum(x * p.reshape(len(p), 1), axis=0)

    def _residuals(self, p, x, y):
        """Return residual between sum of the trends and a light curve."""
        return y - self._func_trends(p, x)

    def detrend(self, lc):
        """Detrend a light curves using the constructed master trends."""

        # Convert the light curve set to numpy array.
        if type(lc) != np.ndarray:
            lc = np.array(lc)

        # Normalize.
        med_lc = np.median(lc)
        raw = (lc - med_lc) / med_lc

        # Initial guess.
        p0 = np.ones(len(self.master_trends))
        # Bounds in [0, infinite]
        p1 = least_squares(self._residuals, p0, args=(self.master_trends, raw),
                           bounds=[0, np.inf])
        p1 = p1['x']

        detrended = raw - self._func_trends(p1, self.master_trends)

        # Scale back to the original flux.
        detrended = detrended * med_lc + med_lc

        return detrended

    def plot_spatial(self, filename='spatial.png'):
        """
        Plot a spatial distribution of the constructed master trends.

        Parameters
        ----------
        filename : str, optional
            A png filename including the path. For example,
            "./outputs/spatial.png". Default is "spatial.png"
        """
        if self.xy_coords is None:
            raise RuntimeError('No x and y coordinates are given.')

        import pylab as pl

        pl.figure(figsize=(12, 12))
        pl.title('Spatial distribution of the constructed master trends')

        colors = 'bgrkmc'
        marks = 's^.+x*'
        for i in range(len(self.master_trends_indices)):
            indices = self.master_trends_indices[i]
            pl.plot(self.xy_coords[indices][:, 0],
                    self.xy_coords[indices][:, 1],
                    marker=marks[int(i / len(colors)) % len(marks)],
                    color=colors[i % len(colors)], ls='None',
                    label='Master trend %d: %d light curves' %
                          (i + 1, len(self.master_trends_indices[i])))
        pl.xlabel('X coordinate')
        pl.ylabel('Y coordinate')
        pl.legend(numpoints=1)
        pl.savefig(filename)
