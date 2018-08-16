"""Tests for `pyssa` package."""

import numpy as np
import pytest

from pyssa.pyssa import direct_naive

def test_null():
    V_r = np.array([[1,0,0],[0,1,0]])
    V_p = np.array([[0,1,0],[0,0,1]])
    X0 = np.array([100,0,0])
    k = np.array([0,0])
    print(V_r, X0, V_r.shape, X0.shape)
    with pytest.raises(RuntimeWarning):
        direct_naive(V_r, V_p, X0, k, max_t = 10, max_iter = 100)

def test_higher_order():
    V_r = np.array([[2,0,0],[0,1,0]])
    V_p = np.array([[0,1,0],[0,0,1]])
    X0 = np.array([100,0,0])
    k = np.array([1,1])
    print(V_r, X0, V_r.shape, X0.shape)
    with pytest.raises(RuntimeWarning):
        direct_naive(V_r, V_p, X0, k, max_t = 10, max_iter = 100)

def test_too_high_order():
    V_r = np.array([[2,2,0],[0,1,0]])
    V_p = np.array([[0,1,0],[0,0,1]])
    X0 = np.array([100,0,0])
    k = np.array([1,1])
    print(V_r, X0, V_r.shape, X0.shape)
    with pytest.raises(RuntimeError):
        direct_naive(V_r, V_p, X0, k, max_t = 10, max_iter = 100)

def test_status_3():
    V_r = np.array([[1,0,0],[0,1,0]])
    V_p = np.array([[0,0,0],[0,0,1]])
    X0 = np.array([10,0,0])
    k = np.array([1,1])
    print(V_r, X0, V_r.shape, X0.shape)
    [_, _, status] = direct_naive(V_r, V_p, X0, k, max_t = 10, max_iter = 100)
    assert status == 3

def test_status_2():
    V_r = np.array([[1,0,0],[0,1,0]])
    V_p = np.array([[0,1,0],[0,0,1]])
    X0 = np.array([10,0,0])
    k = np.array([1,1])
    print(V_r, X0, V_r.shape, X0.shape)
    [_, _, status] = direct_naive(V_r, V_p, X0, k, max_t = 1, max_iter = 100)
    assert status == 2
