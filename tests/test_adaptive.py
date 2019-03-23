"""
    Tests for tau leaping with adaptive step size selection function
"""

import numpy as np
import pytest

from pyssa.tau_adaptive import get_HOR, tau_adaptive, step1, step2, step5
from pyssa.utils import get_kstoc
from pyssa.simulation import Simulation


def test_HOR():
    # A ->
    react_stoic = np.array([[1, 0]])
    assert get_HOR(react_stoic) == 1
    # A ->
    # A ->
    # A ->
    react_stoic = np.array([[1, 1, 1], [0, 0, 0]])
    assert np.all(get_HOR(react_stoic) == [1, 0])
    # A + B ->
    # A ->
    react_stoic = np.array([[1, 1], [1, 0]])
    assert np.all(get_HOR(react_stoic) == [2, 2])
    # A + B ->
    # C ->
    react_stoic = np.array([[1, 0], [1, 0], [0, 1]])
    assert np.all(get_HOR(react_stoic) == [2, 2, 1])
    # A + A ->
    # A ->
    # B ->
    react_stoic = np.array([[2, 1, 0], [0, 0, 1]])
    assert np.all(get_HOR(react_stoic) == [-2, 1])
    # A + A + B ->
    # A + B ->
    # B ->
    react_stoic = np.array([[2, 1, 0], [1, 1, 1]])
    assert np.all(get_HOR(react_stoic) == [-32, 3])
    # A + A + B ->
    # A + B + B ->
    react_stoic = np.array([[2, 1], [1, 2]])
    assert np.all(get_HOR(react_stoic) == [-32, -32])
    # A + A + A ->
    # A + A + B ->
    # A + B ->
    react_stoic = np.array([[3, 2, 1], [0, 1, 1]])
    assert np.all(get_HOR(react_stoic) == [-3, 3])
    # A + A + A ->
    # B + B + B ->
    # A + A + B ->
    # A + B + B ->
    # A + B + C ->
    react_stoic = np.array([[3, 0, 2, 1, 1], [0, 3, 1, 2, 1], [0, 0, 0, 0, 1]])
    assert np.all(get_HOR(react_stoic) == [-3, -3, 3])
    # A + B ->
    # A + A + A ->
    # B + B ->
    # A + C ->
    # C + D + E ->
    react_stoic = np.array(
        [
            [1, 3, 0, 1, 0],
            [1, 0, 2, 0, 0],
            [0, 0, 0, 1, 1],
            [0, 0, 0, 0, 1],
            [0, 0, 0, 0, 1],
        ]
    )
    assert np.all(get_HOR(react_stoic) == [-3, -2, 3, 3, 3])


def test_len_sim(setup_basic):
    V_r, V_p, X0, k = setup_basic
    X0 = X0 = np.array([50, 2, 0])
    sim = Simulation(V_r, V_p, X0, k)
    sim.simulate(algorithm="tau_adaptive")
    x = sim.results.x_list
    t = sim.results.t_list
    assert ~np.all(np.array(x[0][0, :]) == np.array(x[0][1, :]))
    assert ~np.all(np.array(t[0][0]) == np.array(t[0][1]))


def test_tauadaptive_step1(setup_basic):
    V_r, V_p, X0, k = setup_basic
    X0 = np.array([100, 9, 0])
    V = V_p - V_r
    kstoc = get_kstoc(V_r, k, 1.0, False)
    prop, crit, not_crit = step1(kstoc, X0, V_r, V, nc=10)
    assert (prop == [X0[0], X0[1]]).all()
    assert (crit != not_crit).all()
    assert (crit == np.array([False, True])).all()


def test_tauadaptive_step2(setup_basic):
    V_r, V_p, X0, k = setup_basic
    V = V_p - V_r
    HOR = get_HOR(V_r)
    react_species = np.where(np.sum(V_r, axis=1) > 0)[0]
    assert (react_species == np.array([0, 1])).all()
    kstoc = get_kstoc(V_r, k, 1.0, False)
    prop, _, not_crit = step1(kstoc, X0, V_r, V, nc=10)
    taup = step2(not_crit, react_species, V, X0, HOR, prop, epsilon=0.03)
    assert taup == 0.01


def test_tauadaptive_step5(setup_basic):
    V_r, V_p, X0, k = setup_basic
    nr = V_r.shape[1]
    V = V_p - V_r
    kstoc = get_kstoc(V_r, k, 1.0, False)
    prop, _, not_crit = step1(kstoc, X0, V_r, V, nc=10)
    taup, taupp = 100, 10
    tau, K = step5(taup, taupp, nr, not_crit, prop, X0)
    assert tau == taupp
    assert K[1] == 0
    taup, taupp = 0.01, 1
    tau, K = step5(taup, taupp, nr, not_crit, prop, X0)
    assert tau == taup
    assert K[1] == 0
