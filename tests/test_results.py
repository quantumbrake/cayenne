"""
    Tests for the `Results` class
"""

import pytest

from pyssa.simulation import Simulation
from pyssa.results import Results


@pytest.mark.parametrize("algorithm", ["direct", "tau_leaping"])
@pytest.mark.usefixtures("setup_large")
class TestResults:
    def test_init_good(self, algorithm, setup_large):
        V_r, V_p, X0, k = setup_large
        sim = Simulation(V_r, V_p, X0, k)
        with pytest.warns(Warning):
            sim.results
        sim.simulate(algorithm=algorithm)
        assert sim.results

    def test_init_bad(self, algorithm, setup_large):
        V_r, V_p, X0, k = setup_large
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm)
        results = sim.results
        t_list = results.t_list
        x_list = results.x_list
        status_list = results.status_list
        algorithm = "direct"
        seed = [0] * len(status_list)
        assert Results(t_list, x_list, status_list, algorithm, seed)
        with pytest.raises(ValueError):
            Results(t_list, x_list, status_list, algorithm, seed[:-2])
        with pytest.raises(ValueError):
            Results(t_list[:-2], x_list, status_list, algorithm, seed)
        with pytest.raises(ValueError):
            t_list[0] = t_list[0][:-2]
            Results(t_list, x_list, status_list, algorithm, seed)
        with pytest.raises(ValueError):
            status_list[-1] = "fail"
            Results(t_list, x_list, status_list, algorithm, seed)

    def test_iter_len(self, algorithm, setup_large):
        V_r, V_p, X0, k = setup_large
        n_rep = 10
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, n_rep=n_rep)
        results = sim.results
        assert len(results) == n_rep
        for x, t, s in results:
            assert x.shape[0] == t.shape[0]
            assert s

    def test_contains_getitem(self, algorithm, setup_large):
        V_r, V_p, X0, k = setup_large
        n_rep = 10
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, n_rep=n_rep)
        results = sim.results
        assert 9 in results
        x, t, s = results[9]
        assert x.shape[0] == t.shape[0]
        assert s
        with pytest.raises(IndexError):
            results[100]

    def test_final(self, algorithm, setup_large):
        V_r, V_p, X0, k = setup_large
        n_rep = 3
        max_t = 1e5
        max_iter = 100
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, max_t=max_t, max_iter=max_iter, n_rep=n_rep)
        results = sim.results
        final_times, final_states = results.final
        assert final_times.shape[0] == final_states.shape[0]
        if algorithm != "tau_leaping":
            for i in range(final_states.shape[0]):
                assert (final_states[0, :] == final_states[i, :]).all()
