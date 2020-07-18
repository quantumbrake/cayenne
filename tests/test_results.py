"""
    Tests for the `Results` class
"""

import numpy as np
import pytest

from cayenne.simulation import Simulation
from cayenne.results import Results


@pytest.mark.parametrize("algorithm", ["direct", "tau_leaping", "tau_adaptive"])
@pytest.mark.usefixtures("setup_large")
class TestResults:
    """
        Test the Results class.
    """

    def test_init_good(self, algorithm, setup_large):
        """
            Test if initialization works.
        """
        species_names, rxn_names, V_r, V_p, X0, k = setup_large
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        with pytest.warns(Warning):
            sim.results
        sim.simulate(algorithm=algorithm)
        assert sim.results

    def test_init_bad(self, algorithm, setup_large):
        """
            Test if initialization fails when it is supposed to.
        """
        species_names, rxn_names, V_r, V_p, X0, k = setup_large
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm)
        results = sim.results
        t_list = results.t_list
        x_list = results.x_list
        status_list = results.status_list
        algorithm = "direct"
        seed = [0] * len(status_list)
        assert Results(
            species_names, rxn_names, t_list, x_list, status_list, algorithm, seed
        )
        with pytest.raises(ValueError):
            Results(
                species_names,
                rxn_names,
                t_list,
                x_list,
                status_list,
                algorithm,
                seed[:-2],
            )
        with pytest.raises(ValueError):
            Results(
                species_names,
                rxn_names,
                t_list[:-2],
                x_list,
                status_list,
                algorithm,
                seed,
            )
        with pytest.raises(ValueError):
            t_list[0] = t_list[0][:-2]
            Results(
                species_names, rxn_names, t_list, x_list, status_list, algorithm, seed
            )
        with pytest.raises(ValueError):
            status_list[-1] = "fail"
            Results(
                species_names, rxn_names, t_list, x_list, status_list, algorithm, seed
            )

    def test_iter_len(self, algorithm, setup_large):
        """
            Test the ``__len__`` method and if shape of x and t match.
        """
        species_names, rxn_names, V_r, V_p, X0, k = setup_large
        n_rep = 10
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, n_rep=n_rep)
        results = sim.results
        assert len(results) == n_rep
        for x, t, s in results:
            assert x.shape[0] == t.shape[0]
            assert s

    def test_contains_getitem(self, algorithm, setup_large):
        """
            Test the ``__getitem__`` method.
        """
        species_names, rxn_names, V_r, V_p, X0, k = setup_large
        n_rep = 10
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, n_rep=n_rep)
        results = sim.results
        assert 9 in results
        x, t, s = results[9]
        assert x.shape[0] == t.shape[0]
        assert s
        with pytest.raises(IndexError):
            results[100]

    def test_final(self, algorithm, setup_large):
        """
            Test the ``final`` property.
        """
        species_names, rxn_names, V_r, V_p, X0, k = setup_large
        n_rep = 3
        max_t = 1e5
        max_iter = 100
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, max_t=max_t, max_iter=max_iter, n_rep=n_rep)
        results = sim.results
        final_times, final_states = results.final
        assert final_times.shape[0] == final_states.shape[0]
        if algorithm != "tau_leaping" and algorithm != "tau_leaping":
            for i in range(final_states.shape[0]):
                assert (final_states[0, :] == final_states[i, :]).all()

    def test_get_states(self, algorithm, setup_large):
        """
            Test the ``get_state`` method.
        """
        species_names, rxn_names, V_r, V_p, X0, k = setup_large
        n_rep = 3
        max_t = 1e5
        max_iter = 100
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, max_t=max_t, max_iter=max_iter, n_rep=n_rep)
        results = sim.results
        assert np.isclose(results.get_state(0.0), np.array(X0, dtype=np.float)).all()

        x_list = [
            np.array([0, 1, 2, 3]).reshape(-1, 1),
            np.array([0, 2, 3, 4]).reshape(-1, 1),
            np.array([0, 0, 1, 2]).reshape(-1, 1),
        ]
        t_list = [
            np.array([0, 1, 2, 3]),
            np.array([0, 1, 2, 3]),
            np.array([0, 1, 2, 3]),
        ]
        status_list = [1, 1, 1]
        seed = [0, 1, 2]
        res = Results(["A"], ["r1"], t_list, x_list, status_list, algorithm, seed)
        zero_array = np.array([0])
        one_array = np.array([1])
        two_array = np.array([2])
        three_array = np.array([3])
        four_array = np.array([4])
        assert np.isclose(res.get_state(0), [zero_array, zero_array, zero_array]).all()
        assert np.isclose(res.get_state(1), [one_array, two_array, zero_array]).all()
        assert np.isclose(
            res.get_state(1 + np.finfo(float).eps), [one_array, two_array, zero_array]
        ).all()
        assert np.isclose(res.get_state(3), [three_array, four_array, two_array]).all()

    def test_get_states_high_time(self, algorithm, setup_00003):
        """
            Test the ``get_state`` method at longer model times.
        """
        (
            species_names,
            rxn_names,
            V_r,
            V_p,
            X0,
            k,
            _,
            _,
            _,
            max_t,
            max_iter,
            n_rep,
        ) = setup_00003
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, max_t=10000, max_iter=max_iter, n_rep=1)
        print(
            "status", sim.results.status_list[0], sim.results.t_list, sim.results.x_list
        )
        states = sim.results.get_state(t=15000)
        _, x_final = sim.results.final
        assert states[0] == x_final[0]
        sim.simulate(algorithm=algorithm, max_t=2, max_iter=max_iter, n_rep=1)
        with pytest.warns(UserWarning):
            sim.results.get_state(t=1_000_000)
