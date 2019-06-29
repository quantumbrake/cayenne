"""Tests for the cython implementation of direct_naive."""

import pytest
import numpy as np
from pyssa.simulation import Simulation
from pyssa.utils_cython import roulette_selection, get_kstoc


@pytest.mark.usefixtures("setup_basic", "setup_large")
class TestCython:
    def test_sum(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(
            algorithm="direct_cython", max_t=1e5, max_iter=int(1e8), chem_flag=False
        )

    def test_get_kstoc(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        kstoc = get_kstoc(V_r, k, 1.0, False)
        assert (k == kstoc).all()
        kstoc = get_kstoc(V_r, k, 2.0, False)
        assert (k == kstoc).all()
        kstoc = get_kstoc(V_r, k, 2.0, True)
        assert (k == kstoc).all()


class TestRoulette:
    def test_roulette(self):
        choice, status = roulette_selection(np.array([5, 0, 0]), [1, 2])
        assert choice == 0
        assert status == 0
        choice, status = roulette_selection(np.array([0, 5, 0]), [1, 2])
        assert choice == 1
        assert status == 0
        choice, status = roulette_selection(np.array([0, 0, 5]), [1, 2])
        assert choice == 2
        assert status == 0
        choice, status = roulette_selection(np.array([0, 0, 0]), [0, 0])
        assert choice == -1
        assert status == 3
        choice, status = roulette_selection(np.array([0, 0, 0]), [1, 2])
        assert choice == -1
        assert status == -2
