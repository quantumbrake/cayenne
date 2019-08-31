"""
    Tests for tau leaping with adaptive step size selection function
"""

import numpy as np
from pyssa.algorithms.tau_adaptive_cython import py_step1
from pyssa.utils_cython import get_kstoc


def test_tauadaptive_step1(setup_basic):
    V_r, V_p, X0, k = setup_basic
    V = V_p - V_r
    kstoc = get_kstoc(V_r, k, 1.0, False)
    X0 = np.array([100, 9, 0], dtype=np.int64)
    prop, crit, not_crit = py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit != not_crit).all()
    assert (crit == np.array([False, True])).all()

    X0 = np.array([100, 11, 0], dtype=np.int64)
    prop, crit, not_crit = py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit == np.array([False, False])).all()

    X0 = np.array([100, 0, 0], dtype=np.int64)
    prop, crit, not_crit = py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit == np.array([False, False])).all()

    X0 = np.array([0, 9, 0], dtype=np.int64)
    prop, crit, not_crit = py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit == np.array([False, True])).all()

    X0 = np.array([0, 0, 0], dtype=np.int64)
    prop, crit, not_crit = py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit == np.array([False, False])).all()

    # 0 -> A
    # A -> 0
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    V = V_p - V_r
    X0 = np.array([100], dtype=np.int64)
    prop, crit, not_crit = py_step1(kstoc, X0, V_r, V, nc=10)
    assert (crit == np.array([False, False])).all()

    X0 = np.array([9], dtype=np.int64)
    prop, crit, not_crit = py_step1(kstoc, X0, V_r, V, nc=10)
    assert (crit == np.array([False, True])).all()

    X0 = np.array([0], dtype=np.int64)
    prop, crit, not_crit = py_step1(kstoc, X0, V_r, V, nc=10)
    assert (crit == np.array([False, False])).all()


# def test_len_sim(setup_basic):
#     V_r, V_p, X0, k = setup_basic
#     X0 = X0 = np.array([50, 2, 0])
#     sim = Simulation(V_r, V_p, X0, k)
#     sim.simulate(algorithm="tau_adaptive")
#     x = sim.results.x_list
#     t = sim.results.t_list
#     assert ~np.all(np.array(x[0][0, :]) == np.array(x[0][1, :]))
#     assert ~np.all(np.array(t[0][0]) == np.array(t[0][1]))


# def test_tauadaptive_step2(setup_basic):
#     V_r, V_p, X0, k = setup_basic
#     V = V_p - V_r
#     HOR = get_HOR(V_r)
#     react_species = np.where(np.sum(V_r, axis=1) > 0)[0]
#     assert (react_species == np.array([0, 1])).all()
#     kstoc = get_kstoc(V_r, k, 1.0, False)
#     prop, _, not_crit = step1(kstoc, X0, V_r, V, nc=10)
#     taup = step2(not_crit, react_species, V, X0, HOR, prop, epsilon=0.03)
#     assert taup == 0.01


# def test_tauadaptive_step5(setup_basic):
#     V_r, V_p, X0, k = setup_basic
#     nr = V_r.shape[1]
#     V = V_p - V_r
#     kstoc = get_kstoc(V_r, k, 1.0, False)
#     prop, _, not_crit = step1(kstoc, X0, V_r, V, nc=10)
#     taup, taupp = 100, 10
#     tau, K = step5(taup, taupp, nr, not_crit, prop, X0)
#     assert tau == taupp
#     assert K[1] == 0
#     taup, taupp = 0.01, 1
#     tau, K = step5(taup, taupp, nr, not_crit, prop, X0)
#     assert tau == taup
#     assert K[1] == 0
