"""
    Tests for kstoc function
"""

import numpy as np
import pytest

from pyssa.utils import Na, get_kstoc
from pyssa.simulation import Simulation


@pytest.mark.usefixtures("setup_system")
class TestKstoc:
    def test_100(self, setup_system):
        k_det, volume = setup_system
        V_r = np.array([[1], [0], [0]])
        k_stoc = get_kstoc(V_r, k_det, volume, chem_flag=True)
        assert k_stoc == k_det

    def test_110(self, setup_system):
        k_det, volume = setup_system
        V_r = np.array([[1], [1], [0]])
        k_stoc = get_kstoc(V_r, k_det, volume, chem_flag=True)
        assert k_stoc == k_det / (Na * volume)

    def test_111(self, setup_system):
        k_det, volume = setup_system
        V_r = np.array([[1], [1], [1]])
        k_stoc = get_kstoc(V_r, k_det, volume, chem_flag=True)
        assert k_stoc == k_det / (Na * volume) ** 2

    def test_200(self, setup_system):
        k_det, volume = setup_system
        V_r = np.array([[2], [0], [0]])
        k_stoc = get_kstoc(V_r, k_det, volume, chem_flag=True)
        assert k_stoc == k_det * 2 / (Na * volume)

    def test_210(self, setup_system):
        k_det, volume = setup_system
        V_r = np.array([[2], [1], [0]])
        k_stoc = get_kstoc(V_r, k_det, volume, chem_flag=True)
        assert k_stoc == k_det * 2 / (Na * volume) ** 2

    def test_300(self, setup_system):
        k_det, volume = setup_system
        V_r = np.array([[3], [0], [0]])
        k_stoc = get_kstoc(V_r, k_det, volume, chem_flag=True)
        assert k_stoc == k_det * 6 / (Na * volume) ** 2

    def test_chemflag(self, setup_system):
        k_det, volume = setup_system
        V_r = np.array([[3], [0], [0]])
        k_stoc = get_kstoc(V_r, k_det, volume, chem_flag=True)
        assert k_stoc == k_det * 6 / (Na * volume) ** 2
        k_det, volume = setup_system
        V_r = np.array([[3], [0], [0]])
        k_stoc = get_kstoc(V_r, k_det, volume, chem_flag=False)
        assert k_stoc == k_det * 6 / (1 * volume) ** 2
