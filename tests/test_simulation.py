"""
    Tests for direct_naive function
"""

import numpy as np
import pytest

from pyssa.simulation import Simulation


@pytest.mark.parametrize("algorithm", ["direct_naive", "tau_leaping"])
@pytest.mark.usefixtures("setup_basic", "setup_large")
class TestSanitize:
    def test_null(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        k = np.array([0, 0])
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm)
        assert sim.results.status_list[0] == -2

    def test_too_high_order(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_r = np.array([[2, 2, 0], [0, 1, 0]])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_status_3(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, 0, 0], [0, 0, 1]])
        X0 = np.array([10, 0, 0])
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, max_t=10, max_iter=100, chem_flag=True)
        assert sim.results.status_list[0] == 3

    def test_status_2(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        X0 = np.array([10, 0, 0])
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, max_t=1, max_iter=100, chem_flag=True)
        assert sim.results.status_list[0] == 2

    def test_neg_k(self, algorithm, setup_large):
        V_r, V_p, X0, k = setup_large
        k = np.array([1, 1, -1, 1, -1])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_Vp_Vr_shape(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, 1, 0]])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_kdet_Vr_shape(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        k = np.array([1, 1, 1])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_Vp_neg(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, -1, 0], [0, 0, 1]])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_Vr_neg(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_r = np.array([[-1, 0, 0], [0, 1, 0]])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_X0_neg(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        X0 = np.array([-10, 0, 0])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_reproduce(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(V_r, V_p, X0, k)
        sim2 = Simulation(V_r, V_p, X0, k)
        sim1.simulate(algorithm=algorithm)
        sim2.simulate(algorithm=algorithm)
        assert all(
            (i == j).all() for i, j in zip(sim1.results.t_list, sim2.results.t_list)
        )
        assert all(
            (i == j).all() for i, j in zip(sim1.results.x_list, sim2.results.x_list)
        )
        assert all(
            i == j for i, j in zip(sim1.results.status_list, sim2.results.status_list)
        )

    def test_reproduce_fail(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(V_r, V_p, X0, k)
        sim2 = Simulation(V_r, V_p, X0, k)
        sim1.simulate(algorithm=algorithm)
        sim2.simulate(algorithm=algorithm, seed=[1])
        for i, j in zip(sim1.results.t_list, sim2.results.t_list):
            if i.shape[0] == j.shape[0]:
                assert not (
                    all(
                        (i == j).all()
                        for i, j in zip(sim1.results.t_list, sim2.results.t_list)
                    )
                )

    def test_incorrect_seed(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(V_r, V_p, X0, k)
        with pytest.raises(ValueError):
            sim1.simulate(algorithm=algorithm, n_rep=2, seed=[1])
        with pytest.raises(ValueError):
            sim1.simulate(algorithm=algorithm, n_rep=2, seed=[1, 2, 3])

    def test_maxiter_type(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(V_r, V_p, X0, k)
        with pytest.raises(TypeError):
            sim1.simulate(algorithm=algorithm, max_iter=100.0)

    def test_tau(self, algorithm, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(algorithm="tau_leaping", n_rep=1)
        results = sim.results
        for x, t, s in results:
            print(s, x)
        assert 1 == 1


def test_bifurcation(setup_bifurcation):
    V_r, V_p, X0, k = setup_bifurcation
    count_excitation = 0
    n_runs = 10
    deviation_tolerance = 0.05
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=n_runs)
    for ind in range(n_runs):
        x = sim1.results.x_list[ind]
        xt = x[-1, :]
        assert np.all(xt - np.array([0, 11, 0, 0]) == 0) or np.all(
            xt - np.array([0, 0, 1, 10]) == 0
        )
        if np.all(xt - np.array([0, 11, 0, 0]) == 0):
            count_excitation += 1
    assert np.abs(count_excitation / n_runs - 0.5) < deviation_tolerance


# def test_long(setup_long):
#     V_r, V_p, X0, k = setup_long
#     _, Xt, status = direct_naive(
#         V_r, V_p, X0, k, max_t=1e5, max_iter=1e8, chem_flag=False
#     )
#     X_output = np.array([0, 0, X0[0]])
#     assert status == -2
#     assert Xt.all() == X_output.all()
