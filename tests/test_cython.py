"""Tests for the cython implementation of direct_naive."""

import pytest
import numpy as np
from cayenne.simulation import Simulation
from cayenne.utils import get_kstoc
from cayenne.utils import py_roulette_selection as roulette_selection


@pytest.mark.usefixtures("setup_basic", "setup_large")
class TestCython:
    def test_sim(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm="direct", max_t=1e5, max_iter=int(1e8), chem_flag=False)

    def test_get_kstoc(self, setup_basic):

        # A -> B
        # B -> C
        _, _, vr, _, _, k = setup_basic
        kstoc = np.array(get_kstoc(vr, k, 1.0, False))
        assert (kstoc == k).all()
        kstoc = np.array(get_kstoc(vr, k, 2.0, False))
        assert (kstoc == k).all()
        kstoc = np.array(get_kstoc(vr, k, 2.0, True))
        assert (kstoc == k).all()

        # A + B -> C
        # B + C -> C
        vr = np.array([[1, 1, 0], [1, 1, 0]]).transpose()
        kstoc = np.array(get_kstoc(vr, k, 1.0, False))
        V = 3.0
        print(kstoc)
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = np.array(get_kstoc(vr, k, V, False))
        assert (kstoc == [k[0] / V, k[1] / V]).all()

        # A + B -> C
        # A -> C
        vr = np.array([[1, 1, 0], [1, 0, 0]]).transpose()
        kstoc = np.array(get_kstoc(vr, k, 1.0, False))
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = np.array(get_kstoc(vr, k, V, False))
        assert (kstoc == [k[0] / V, k[1]]).all()

        # A -> C
        # A + B -> C
        vr = np.array([[1, 0, 0], [1, 1, 0]]).transpose()
        kstoc = np.array(get_kstoc(vr, k, 1.0, False))
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = np.array(get_kstoc(vr, k, V, False))
        assert (kstoc == [k[0], k[1] / V]).all()

        # A + B -> C
        vr = np.array([[1, 1, 0]]).transpose()
        k = np.array([1.0])
        kstoc = np.array(get_kstoc(vr, k, 1.0, False))
        assert kstoc == k
        kstoc = np.array(get_kstoc(vr, k, V, False))
        assert kstoc == k / V

        # A + B + C -> D
        # A -> B
        vr = np.array([[1, 1, 1, 0], [1, 0, 0, 0]]).transpose()
        k = np.array([1.0, 1.0])
        kstoc = np.array(get_kstoc(vr, k, 1.0, False))
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = np.array(get_kstoc(vr, k, V, False))
        assert (kstoc == [k[0] / V ** 2, k[1]]).all()

        # A + B -> C
        # A + B + C -> D
        vr = np.array([[1, 1, 0, 0], [1, 1, 1, 0]]).transpose()
        k = np.array([1.0, 1.0])
        kstoc = np.array(get_kstoc(vr, k, 1.0, False))
        assert (kstoc == [k[0], k[1]]).all()
        kstoc = np.array(get_kstoc(vr, k, V, False))
        assert (kstoc == [k[0] / V, k[1] / V ** 2]).all()


class TestRoulette:
    def test_roulettecy(self):
        choice, status = roulette_selection(
            np.array([5, 0, 0], dtype=np.double), np.array([1, 2], dtype=np.int64)
        )
        assert choice == 0
        assert status == 0
        choice, status = roulette_selection(
            np.array([0, 5, 0], dtype=np.double), np.array([1, 2], dtype=np.int64)
        )
        assert choice == 1
        assert status == 0
        choice, status = roulette_selection(
            np.array([0, 0, 5], dtype=np.double), np.array([1, 2], dtype=np.int64)
        )
        assert choice == 2
        assert status == 0
        choice, status = roulette_selection(
            np.array([0, 0, 0], dtype=np.double), np.array([0, 0], dtype=np.int64)
        )
        assert choice == -1
        assert status == 3
        choice, status = roulette_selection(
            np.array([0, 0, 0], dtype=np.double), np.array([1, 2], dtype=np.int64)
        )
        assert choice == -1
        assert status == -2

    def test_100(self):
        prop_list = np.array([1.0, 0.0, 0.0])
        xt = np.array([1, 2], dtype=np.int64)
        choice, status = roulette_selection(prop_list, xt)
        assert status == 0
        assert choice == 0

    def test_010(self):
        prop_list = np.array([0.0, 1.0, 0.0])
        xt = np.array([1, 2], dtype=np.int64)
        choice, status = roulette_selection(prop_list, xt)
        assert status == 0
        assert choice == 1

    def test_001(self):
        prop_list = np.array([0.0, 0.0, 1.0])
        xt = np.array([1, 2], dtype=np.int64)
        choice, status = roulette_selection(prop_list, xt)
        assert status == 0
        assert choice == 2

    def test_stat3(self):
        prop_list = np.array([0.0, 0.0, 0.0])
        xt = np.array([0, 0], dtype=np.int64)
        _, status = roulette_selection(prop_list, xt)
        assert status == 3

    def test_statm2(self):
        prop_list = np.array([0.0, 0.0, 0.0])
        xt = np.array([1, 0], dtype=np.int64)
        _, status = roulette_selection(prop_list, xt)
        assert status == -2
