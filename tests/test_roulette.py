"""
    Tests for roulette_selection function
"""

import numpy as np

from pyssa.utils import roulette_selection


class TestRoulette:
    def test_100(self):
        prop_list = np.array([1, 0, 0])
        xt = np.array([1, 2])
        choice, status = roulette_selection(prop_list, xt)
        assert status == 0
        assert choice == 0

    def test_010(self):
        prop_list = np.array([0, 1, 0])
        xt = np.array([1, 2])
        choice, status = roulette_selection(prop_list, xt)
        assert status == 0
        assert choice == 1

    def test_001(self):
        prop_list = np.array([0, 0, 1])
        xt = np.array([1, 2])
        choice, status = roulette_selection(prop_list, xt)
        assert status == 0
        assert choice == 2

    def test_stat3(self):
        prop_list = np.array([0, 0, 0])
        x_t = np.array([0, 0])
        _, status = roulette_selection(prop_list, x_t)
        assert status == 3

    def test_statm2(self):
        prop_list = np.array([0, 0, 0])
        x_t = np.array([1, 0])
        _, status = roulette_selection(prop_list, x_t)
        assert status == -2
