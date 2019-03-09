"""
    Tests for tau leaping with adaptive step size selection function
"""

import numpy as np
import pytest

from pyssa.tau_adaptive import get_HOR


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
