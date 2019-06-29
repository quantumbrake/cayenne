"""Tests for the cython implementation of direct_naive."""

import pytest
import numpy as np
from pyssa.simulation import Simulation
from pyssa.utils_cython import roulette_selection, get_kstoc


@pytest.mark.usefixtures("setup_basic", "setup_large")
class TestCython:
    def test_sim(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(
            algorithm="direct_cython", max_t=1e5, max_iter=int(1e8), chem_flag=False
        )

    def test_get_kstoc(self, setup_basic):

        # A -> B
        # B -> C
        vr, _, _, k = setup_basic
        kstoc = get_kstoc(vr, k, 1.0, False)
        assert (kstoc == k).all()
        kstoc = get_kstoc(vr, k, 2.0, False)
        assert (kstoc == k).all()
        kstoc = get_kstoc(vr, k, 2.0, True)
        assert (kstoc == k).all()

        # A + B -> C
        # B + C -> C
        vr = np.array([[1, 1, 0], [1, 1, 0]]).transpose()
        kstoc = get_kstoc(vr, k, 1.0, False)
        V = 3.0
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = get_kstoc(vr, k, V, False)
        assert (kstoc == [k[0] / V, k[1] / V]).all()

        # A + B -> C
        # A -> C
        vr = np.array([[1, 1, 0], [1, 0, 0]]).transpose()
        kstoc = get_kstoc(vr, k, 1.0, False)
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = get_kstoc(vr, k, V, False)
        assert (kstoc == [k[0] / V, k[1]]).all()

        # A -> C
        # A + B -> C
        vr = np.array([[1, 0, 0], [1, 1, 0]]).transpose()
        kstoc = get_kstoc(vr, k, 1.0, False)
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = get_kstoc(vr, k, V, False)
        assert (kstoc == [k[0], k[1] / V]).all()

        # A + B -> C
        vr = np.array([[1, 1, 0]]).transpose()
        k = np.array([1.0])
        kstoc = get_kstoc(vr, k, 1.0, False)
        assert kstoc == k
        kstoc = get_kstoc(vr, k, V, False)
        assert kstoc == k / V

        # A + B + C -> D
        # A -> B
        vr = np.array([[1, 1, 1, 0], [1, 0, 0, 0]]).transpose()
        k = np.array([1.0, 1.0])
        kstoc = get_kstoc(vr, k, 1.0, False)
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = get_kstoc(vr, k, V, False)
        assert (kstoc == [k[0] / V ** 2, k[1]]).all()

        # A + B -> C
        # A + B + C -> D
        vr = np.array([[1, 1, 0, 0], [1, 1, 1, 0]]).transpose()
        k = np.array([1.0, 1.0])
        kstoc = get_kstoc(vr, k, 1.0, False)
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = get_kstoc(vr, k, V, False)
        assert (kstoc == [k[0] / V, k[1] / V ** 2]).all()


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
