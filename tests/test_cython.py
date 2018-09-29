"""
    Tests for cy_direct_naive function
"""

import numpy as np
import pytest

from pyssa.pyssa_cython import cy_direct_naive, Na


@pytest.mark.usefixtures("setup_basic", "setup_large")
class TestSanitize:

    def test_null(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        k = np.array([0.0, 0.0])
        [_, _, status] = cy_direct_naive(V_r, V_p, X0, k, max_t=10, max_iter=100)
        assert status == -2

    def test_too_high_order(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_r = np.array([[2, 2, 0], [0, 1, 0]])
        with pytest.raises(ValueError):
            cy_direct_naive(V_r, V_p, X0, k, max_t=10, max_iter=100)

    def test_status_3(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, 0, 0], [0, 0, 1]])
        X0 = np.array([10, 0, 0])
        [_, _, status] = cy_direct_naive(V_r, V_p, X0, k, max_t=10, max_iter=100)
        assert status == 3

    def test_status_2(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        X0 = np.array([10, 0, 0])
        [_, _, status] = cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)
        assert status == 2

    def test_neg_k(self, setup_large):
        V_r, V_p, X0, k = setup_large
        k = np.array([1, 1, -1, 1, -1])
        with pytest.raises(ValueError):
            cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_Vp_Vr_shape(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, 1, 0]])
        with pytest.raises(ValueError):
            cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_kdet_Vr_shape(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        k = np.array([1, 1, 1])
        with pytest.raises(ValueError):
            cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_Vp_neg(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, -1, 0], [0, 0, 1]])
        with pytest.raises(ValueError):
            cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_Vr_neg(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_r = np.array([[-1, 0, 0], [0, 1, 0]])
        with pytest.raises(ValueError):
            cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_X0_neg(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        X0 = np.array([-10, 0, 0])
        with pytest.raises(ValueError):
            cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100)

    def test_reproduce(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        [t1, Xt1, status1] = cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100, seed=0)
        [t2, Xt2, status2] = cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100, seed=0)
        assert t1 == t2
        assert np.all(Xt1) == np.all(Xt2)
        assert status1 == status2

    def test_reproduce_fail(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        [t1, _, _] = cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100, seed=0)
        [t2, _, _] = cy_direct_naive(V_r, V_p, X0, k, max_t=1, max_iter=100, seed=1)
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
        [_, Xt, _] = cy_direct_naive(V_r, V_p, X0, k, max_t=150, max_iter=1000, seed=ind)
        assert np.all(Xt - np.array([0, 11, 0, 0]) == 0) or np.all(Xt - np.array([0, 0, 1, 10]) == 0)
        if np.all(Xt - np.array([0, 11, 0, 0]) == 0):
            count_excitation += 1
    assert np.abs(count_excitation / n_runs - 0.5) < deviation_tolerance
