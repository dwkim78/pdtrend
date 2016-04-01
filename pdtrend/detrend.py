__author__ = 'kim'

import numpy as np

from sklearn.cluster import Birch
from scipy.optimize import least_squares


class PDTrend:
    def __init__(self, lcs, n_min_member=5, dist_cut=0.45):
        """
        pdtrend determines master trends using Pearson correlation coefficients
        between each light curve and then provide a function to detrend
        each light curve using the determined master trends. For details,
        see Kim et al. 2009.
        :param lcs: A list of light curves. It must be the shape of N x M,
        where N is the number of light curves and M is the number of data points.
        Note that M must be same for all N light curves. pdtrend assumes that
        these M data points are synced in time. Thus pdtrend does not use
        any time information.
        :param n_min_member: The minimum number of members in each cluster.
        :param dist_cut: Distance cut to filter clusters.
        """

        # Convert the light curve set to numpy array.
        if type(lcs) != np.ndarray:
            lcs = np.array(lcs)

        # Set parameters.
        self.lcs = lcs
        self.n_min_member = n_min_member
        self.dist_cut = dist_cut

    def _calculate_distance_matrix(self):
        """
        Calculate a distance matrix, which is defined as:
        (1. - correlation_matrix) / 2.
        """
        corr_matrix = np.corrcoef(self.lcs)
        dist_matrix = (1. - corr_matrix) / 2.

        self.corr_matrix = corr_matrix
        self.dist_matrix = dist_matrix

    def _find_clusters(self, branching_factor=50, threshold=0.1):
        """
        Find clusters using Birch and the distance matrix.
        """
        birch = Birch(branching_factor=branching_factor,
                      threshold=threshold).fit(self.dist_matrix)

        self.birch = birch

    def _filter_clusters(self):
        """
        Filter a cluster 1) which has less than "n_min_member" members,
        2) median distance between each member is larger than "dist_cut".
        """
        unique_labels = set(self.birch.labels_)
        _filtered_labels = []

        for label in unique_labels:
            index = [i for i in range(len(self.birch.labels_)) if
                     self.birch.labels_[i] == label]
            if len(index) < self.n_min_member:
                continue

            dist_list = []
            for i in range(len(index) - 1):
                for j in range(i + 1, len(index)):
                    dist_list.append(self.dist_matrix[index[i], index[j]])

            if np.median(dist_list) <= self.dist_cut:
                _filtered_labels.append(label)

        self._filtered_labels = _filtered_labels

    def _build_master_trends(self):
        """
        Build master trends using the filtered clusters.
        """
        master_trends_indices = []
        master_trends = []
        for label in self._filtered_labels:
            index = [i for i in range(len(self.birch.labels_)) if
                     self.birch.labels_[i] == label]
            master_trends_indices.append(index)

            trends = []
            for i in range(len(index)):
                # Normalize by standard deviation.
                # So, no need to calculate weight again.
                normed = (self.lcs[index[i]] - np.median(self.lcs[index[i]])) / \
                         np.std(self.lcs[index[i]])
                trends.append(normed)

            master_trends.append(np.sum(trends, axis=0) / len(index))

        self.master_trends_indices = master_trends_indices
        self.master_trends = master_trends

    def run(self):
        """
        Run pdtrend pipeline.
        """
        self._calculate_distance_matrix()
        self._find_clusters()
        self._filter_clusters()
        self._build_master_trends()

    def _func_trends(self, p, x):
        """
        Return sum of the trends.
        """
        return np.sum(x * p.reshape(len(p), 1), axis=0)

    def _residuals(self, p, x, y):
        """
        Return residual between sum of the trends and a light curve.
        """
        return y - self._func_trends(p, x)

    def detrend(self, lc):
        """
        Detrend a light curves using the constructed master trends.
        """
        # Normalize.
        raw = (lc - np.median(lc)) / np.std(lc)

        # Initial guess.
        p0 = np.ones(len(self.master_trends))
        # Bounds in [0, infinite]
        p1 = least_squares(self._residuals, p0, args=(self.master_trends, raw),
                           bounds=[0, np.inf])
        p1 = p1['x']

        detrended = raw - self._func_trends(p1, self.master_trends)

        # Scale back to the original flux.
        detrended = detrended * np.std(lc) + np.median(lc)

        return detrended
