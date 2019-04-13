"""
    Accuracy tests of the algorithms
"""
import numpy as np
import pytest

from pyssa.simulation import Simulation


@pytest.mark.parametrize("algorithm", ["direct", "tau_leaping", "tau_adaptive"])
def test_00001(setup_00001, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list = setup_00001
    n_rep = 10
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=51,
        max_iter=int(1.5e3),
        chem_flag=False,
        n_rep=n_rep,
        # debug=True
    )
    assert np.all(sim1.results.final[0] > time_list[-1])
    for ind1, t in enumerate(time_list[1:]):
        results = sim1.results.get_state(t)
        mu_obs = np.mean(results)
        std_obs = np.std(results)
        Z = np.sqrt(n_rep) * (mu_obs - mu_list[ind1 + 1]) / std_list[ind1 + 1]
        Y = np.sqrt(n_rep / 2) * ((std_obs ** 2) / (std_list[ind1 + 1] ** 2) - 1)
        assert -3 < Z < 3
        assert -5 < Y < 5
