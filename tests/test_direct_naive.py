"""Tests for `pyssa` package."""

import numpy as np
import pytest

from pyssa.pyssa import direct_naive, Na


@pytest.fixture
def setup_basic():
    V_r = np.array([[1, 0, 0], [0, 1, 0]])
    V_p = np.array([[0, 1, 0], [0, 0, 1]])
    X0 = np.array([100, 0, 0])
    k = np.array([1, 1])
    return V_r, V_p, X0, k


@pytest.fixture
def setup_large():
    V_r = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
    V_p = np.array([[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [0, 0, 0, 0, 0]])
    X0 = np.array([10, 0, 0, 0, 0])
    k = np.array([1, 1, 1, 1, 1])
    return V_r, V_p, X0, k


class TestSanitize():

    def test_null(self):
        V_r, V_p, X0, k = setup_basic()
        k = np.array([0, 0])
        [_, _, status] = direct_naive(V_r, V_p, X0, k, max_t=10, max_iter=100,
                                      chem_flag=True)
        assert status == -2

    def test_too_high_order(self):
        V_r, V_p, X0, k = setup_basic()
        V_r = np.array([[2, 2, 0], [0, 1, 0]])
        with pytest.raises(ValueError):
            direct_naive(V_r, V_p, X0, k, max_t=10, max_iter=100)

    def test_status_3(self):
        V_r, V_p, X0, k = setup_basic()
        V_p = np.array([[0, 0, 0], [0, 0, 1]])
        X0 = np.array([10, 0, 0])
        [_, _, status] = direct_naive(V_r, V_p, X0, k, max_t=10, max_iter=100,
                                      chem_flag=True)
        assert status == 3

    def test_status_2(self):
        V_r, V_p, X0, k = setup_basic()
        X0 = np.array([10, 0, 0])
        [_, _, status] = direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100,
                                      chem_flag=True)
        assert status == 2

    def test_neg_k(self):
        V_r, V_p, X0, k = setup_large()
        k = np.array([1, 1, -1, 1, -1])
        with pytest.raises(ValueError):
            direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_Vp_Vr_shape(self):
        V_r, V_p, X0, k = setup_basic()
        V_p = np.array([[0, 1, 0]])
        with pytest.raises(ValueError):
            direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_kdet_Vr_shape(self):
        V_r, V_p, X0, k = setup_basic()
        k = np.array([1, 1, 1])
        with pytest.raises(ValueError):
            direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_Vp_neg(self):
        V_r, V_p, X0, k = setup_basic()
        V_p = np.array([[0, -1, 0], [0, 0, 1]])
        with pytest.raises(ValueError):
            direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_Vr_neg(self):
        V_r, V_p, X0, k = setup_basic()
        V_r = np.array([[-1, 0, 0], [0, 1, 0]])
        with pytest.raises(ValueError):
            direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_X0_neg(self):
        V_r, V_p, X0, k = setup_basic()
        X0 = np.array([-10, 0, 0])
        with pytest.raises(ValueError):
            direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_reproduce(self):
        V_r, V_p, X0, k = setup_basic()
        [t1, Xt1, status1] = direct_naive(V_r, V_p, X0, k, max_t=1,
                                          max_iter=100, seed=0, chem_flag=True)
        [t2, Xt2, status2] = direct_naive(V_r, V_p, X0, k, max_t=1,
                                          max_iter=100, seed=0, chem_flag=True)
        assert t1 == t2
        assert Xt1.all() == Xt2.all()
        assert status1 == status2

    def test_reproduce_fail(self):
        V_r, V_p, X0, k = setup_basic()
        [t1, _, _] = direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100,
                                  seed=0, chem_flag=True)
        [t2, _, _] = direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100,
                                  seed=1, chem_flag=True)
        assert t1 != t2


def test_bifurcation():
    V_r = np.array([[1, 0, 0, 0], [0, 1, 0, 1], [1, 0, 0, 0]])
    V_p = np.array([[0, 1, 0, 0], [0, 2, 0, 0], [0, 0, 1, 0]])
    k = np.array([1, 0.01 * Na, 1])
    X0 = np.array([1, 0, 0, 10])
    count_excitation = 0
    n_runs = 1000
    deviation_tolerance = 0.05
    for ind in range(n_runs):
        [_, Xt, _] = direct_naive(V_r, V_p, X0, k, max_t=150, max_iter=1000, seed=ind,
                                  chem_flag=True)
        assert np.all(Xt - np.array([0, 11, 0, 0]) == 0) or np.all(Xt - np.array([0, 0, 1, 10]) == 0)
        if np.all(Xt - np.array([0, 11, 0, 0]) == 0):
            count_excitation += 1
    assert np.abs(count_excitation / n_runs - 0.5) < deviation_tolerance
