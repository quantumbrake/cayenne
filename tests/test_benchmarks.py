"""
    Benchmark tests
"""

from pyssa.pyssa import direct_naive
from pyssa.pyssa_numba import numba_direct_naive


def test_py_benchmark(benchmark, setup_long):
    V_r, V_p, k, X0 = setup_long
    benchmark(direct_naive, V_r, V_p, X0, k, max_t=1e5, max_iter=1e8, seed=0, chem_flag=True)


def test_numba_benchmark(benchmark, setup_long):
    V_r, V_p, k, X0 = setup_long
    benchmark(numba_direct_naive, V_r, V_p, X0, k, max_t=1e5, max_iter=1e8, seed=0, chem_flag=True)
