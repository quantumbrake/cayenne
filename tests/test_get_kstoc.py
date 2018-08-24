"""Tests for `pyssa` package."""

import numpy as np
import pytest

from pyssa.pyssa import get_kstoc, Na

@pytest.fixture
def setup_system():
    k_det = np.array([3.0])
    volume = 7.0
    return k_det, volume

class TestKstoc():

    def test_100(self):
        k_det, volume = setup_system()
        V_r = np.array([[1, 0, 0]])
        k_stoc = get_kstoc(k_det, V_r, volume)
        assert k_stoc == k_det

    def test_110(self):
        k_det, volume = setup_system()
        V_r = np.array([[1, 1, 0]])
        k_stoc = get_kstoc(k_det, V_r, volume)
        assert k_stoc == k_det / (Na * volume)

    def test_111(self):
        k_det, volume = setup_system()
        V_r = np.array([[1, 1, 1]])
        k_stoc = get_kstoc(k_det, V_r, volume)
        assert k_stoc == k_det / (Na * volume) ** 2

    def test_200(self):
        k_det, volume = setup_system()
        V_r = np.array([[2, 0, 0]])
        k_stoc = get_kstoc(k_det, V_r, volume)
        assert k_stoc == k_det * 2 / (Na * volume)

    def test_210(self):
        k_det, volume = setup_system()
        V_r = np.array([[2, 1, 0]])
        k_stoc = get_kstoc(k_det, V_r, volume)
        assert k_stoc == k_det * 2 / (Na * volume) ** 2

    def test_300(self):
        k_det, volume = setup_system()
        V_r = np.array([[3, 0, 0]])
        k_stoc = get_kstoc(k_det, V_r, volume)
        assert k_stoc == k_det * 6 / (Na * volume) ** 2


