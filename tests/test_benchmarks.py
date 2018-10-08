"""
    Benchmark tests
"""

from pyssa.pyssa import direct_naive
from pyssa.pyssa_cython import cy_direct_naive


def test_py_benchmark(benchmark, setup_bifurcation):
    V_r, V_p, k, X0 = setup_bifurcation
    benchmark(direct_naive, V_r, V_p, X0, k, max_t=150, max_iter=1000, seed=0)


def test_cy_benchmark(benchmark, setup_bifurcation):
    V_r, V_p, k, X0 = setup_bifurcation
    benchmark(cy_direct_naive, V_r, V_p, X0, k, max_t=150, max_iter=1000, seed=0)
