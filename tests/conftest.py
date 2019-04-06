"""
    Common configuration for all the tests
"""

import csv
import pathlib

import numpy as np
import pytest

from pyssa.utils import Na


CURR_DIR = pathlib.Path(__file__).parent


@pytest.fixture
def setup_basic():
    V_r = np.array([[1, 0], [0, 1], [0, 0]])
    V_p = np.array([[0, 0], [1, 0], [0, 1]])
    X0 = np.array([100, 0, 0])
    k = np.array([1.0, 1.0])
    return V_r, V_p, X0, k


@pytest.fixture
def setup_large():
    V_r = np.array(
        [
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1],
        ]
    )
    V_p = np.array(
        [
            [0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
        ]
    )
    X0 = np.array([10, 0, 0, 0, 0])
    k = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    return V_r, V_p, X0, k


@pytest.fixture
def setup_system():
    k_det = np.array([3.0])
    volume = 7.0
    return k_det, volume


@pytest.fixture
def setup_bifurcation():
    V_r = np.array([[1, 0, 1], [0, 1, 0], [0, 0, 0], [0, 1, 0]])
    V_p = np.array([[0, 0, 0], [1, 2, 0], [0, 0, 1], [0, 0, 0]])
    k = np.array([1.0, 0.01 * Na, 1.0])
    X0 = np.array([1, 0, 0, 10])
    return V_r, V_p, X0, k


@pytest.fixture
def setup_long():
    V_r = np.array([[1, 0], [0, 1], [0, 0]])
    V_p = np.array([[0, 0], [1, 0], [0, 1]])
    k = np.array([1.0, 1.0])
    X0 = np.array([int(4e5), 1000, 0])
    return V_r, V_p, X0, k


@pytest.fixture
def setup_00001():
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([100])
    k = np.array([0.1, 0.11])
    print(V_r.shape, V_p.shape)
    time_list = []
    mu_list = []
    std_list = []
    with open(f"{CURR_DIR}/data/results_00001.csv") as fid:
        csv_reader = csv.reader(fid, delimiter=",")
        next(csv_reader)
        for time, mu, std in csv_reader:
            time_list.append(float(time))
            mu_list.append(float(mu))
            std_list.append(float(std))
    return V_r, V_p, X0, k, time_list, np.array(mu_list), np.array(std_list)
