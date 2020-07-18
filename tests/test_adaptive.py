"""
    Tests for tau leaping with adaptive step size selection function
"""

import numpy as np
from cayenne.algorithms.tau_adaptive import (
    _py_step1,
    _py_step2,
    _py_step2_get_g,
    _py_step5,
)
from cayenne.utils import get_kstoc, HIGH
from .test_simulation import TestHOR


def test_tauadaptive_step1(setup_basic):
    _, _, V_r, V_p, X0, k = setup_basic
    V = V_p - V_r
    kstoc = get_kstoc(V_r, k, 1.0, False)
    X0 = np.array([100, 9, 0], dtype=np.int64)
    prop, crit, not_crit = _py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit != not_crit).all()
    assert (crit == np.array([False, True])).all()

    X0 = np.array([100, 11, 0], dtype=np.int64)
    prop, crit, not_crit = _py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit == np.array([False, False])).all()

    X0 = np.array([100, 0, 0], dtype=np.int64)
    prop, crit, not_crit = _py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit == np.array([False, False])).all()

    X0 = np.array([0, 9, 0], dtype=np.int64)
    prop, crit, not_crit = _py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit == np.array([False, True])).all()

    X0 = np.array([0, 0, 0], dtype=np.int64)
    prop, crit, not_crit = _py_step1(kstoc, X0, V_r, V, nc=10)
    print(prop, crit, not_crit)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit == np.array([False, False])).all()

    # 0 -> A
    # A -> 0
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    V = V_p - V_r
    X0 = np.array([100], dtype=np.int64)
    prop, crit, not_crit = _py_step1(kstoc, X0, V_r, V, nc=10)
    assert (crit == np.array([False, False])).all()

    X0 = np.array([9], dtype=np.int64)
    prop, crit, not_crit = _py_step1(kstoc, X0, V_r, V, nc=10)
    assert (crit == np.array([False, True])).all()

    X0 = np.array([0], dtype=np.int64)
    prop, crit, not_crit = _py_step1(kstoc, X0, V_r, V, nc=10)
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


def test_tauadaptive_step2(setup_basic):
    _, _, V_r, V_p, X0, k = setup_basic
    V = V_p - V_r
    sim = TestHOR.create_sim_inst(V_r)
    hor = np.int32(sim.HOR)
    assert np.all(hor == [1, 1, 0])
    react_species = np.int32(np.where(np.sum(V_r, axis=1) > 0)[0])

    assert (react_species == np.array([0, 1])).all()
    kstoc = get_kstoc(V_r, k, 1.0, False)
    prop, _, not_crit = _py_step1(kstoc, np.int64(X0), V_r, V, nc=10)
    taup = _py_step2(not_crit, react_species, V, np.int64(X0), hor, prop, epsilon=0.03)
    assert taup == 0.01

    taup = _py_step2(
        np.zeros(not_crit.shape, dtype=np.int32),
        react_species,
        V,
        np.int64(X0),
        hor,
        prop,
        epsilon=0.03,
    )
    assert np.isclose(taup, HIGH)


def test_tauadaptive_step2_get_g():
    assert _py_step2_get_g(1, 10) == 1
    assert _py_step2_get_g(-2, 2) == 3
    assert _py_step2_get_g(-2, 1) == 2
    assert np.isclose(_py_step2_get_g(-2, 10), 2 + (1 / 9))
    assert _py_step2_get_g(-3, 1) == 3.0
    assert np.isclose(_py_step2_get_g(-3, 10), 3 + 1 / 9 + 2 / 8)
    assert _py_step2_get_g(-32, 1) == 3.0
    assert np.isclose(_py_step2_get_g(-32, 10), 3 / 2 * (2 + 1 / 9))


def test_tauadaptive_step5(setup_basic):
    _, _, V_r, V_p, X0, k = setup_basic
    nr = V_r.shape[1]
    V = V_p - V_r
    kstoc = get_kstoc(V_r, k, 1.0, False)
    prop, _, not_crit = _py_step1(kstoc, np.int64(X0), V_r, V, nc=10)
    taup, taupp = 100.0, 10.0
    tau, K = _py_step5(taup, taupp, nr, not_crit, prop, np.int64(X0))
    assert tau == taupp
    assert K[1] == 0
    taup, taupp = 0.01, 1.0
    tau, K = _py_step5(taup, taupp, nr, not_crit, prop, np.int64(X0))
    assert tau == taup
    assert K[1] == 0
