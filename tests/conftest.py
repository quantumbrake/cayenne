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


def read_results(test_id: str):
    filename = f"{CURR_DIR}/data/results_{test_id}.csv"
    time_list = []
    mu_list = []
    std_list = []
    with open(filename) as fid:
        csv_reader = csv.reader(fid, delimiter=",")
        next(csv_reader)
        for time, mu, std in csv_reader:
            time_list.append(float(time))
            mu_list.append(float(mu))
            std_list.append(float(std))
    return time_list, mu_list, std_list


@pytest.fixture
def setup_00001():
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([100])
    k = np.array([0.1, 0.11])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00001")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00003():
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([100])
    k = np.array([1.0, 1.1])
    max_t = 51
    max_iter = int(1e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00003")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00004():
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([10])
    k = np.array([0.1, 0.11])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00004")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00005():
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([10_000])
    k = np.array([0.1, 0.11])
    max_t = 51
    max_iter = int(5e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00005")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00011():
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([100])
    # divide k by 2 because rate expression given in units of concentration
    # in the model file
    k = np.array([0.1, 0.11]) / 2
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00011")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00020():
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    X0 = np.array([0])
    k = np.array([1.0, 0.1])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00020")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00021():
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    X0 = np.array([0])
    k = np.array([10.0, 0.1])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00021")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00022():
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    X0 = np.array([0])
    k = np.array([5.0, 0.1])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00022")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00023():
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    X0 = np.array([0])
    k = np.array([1000.0, 0.1])
    max_t = 51
    max_iter = int(1.5e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00023")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )
