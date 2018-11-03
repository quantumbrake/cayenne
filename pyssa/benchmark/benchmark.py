"""
    Module for benchmarking the performance of various implementations of the Gillespie algorithm
"""

import timeit


def benchmark(setup: str, stmt: str, number: int = 1000000) -> float:
    """
        Function the runs the benchmark

        Parameters
        ----------
        setup : str
            The string containing the initial setup for the profiling
        stmt : str
            The string containing the code to be profiled
        number : int
            The number of iterations to be performed
        *args
            Arguments to be forwarded to the function
        **kwargs
            Keyword arguments to be forwarded to the function

        Returns
        -------
        float
            The time taken to run the `stmt`
    """
    return timeit.timeit(stmt=stmt, setup=setup, number=number)
