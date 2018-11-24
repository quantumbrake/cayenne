"""
    Tests for direct_naive function
"""

import numpy as np
import pytest

from pyssa.simulation import Simulation


@pytest.mark.usefixtures("setup_basic", "setup_large")
class TestSanitize:
    def test_null(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        k = np.array([0, 0])
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate()
        assert sim.results.status_list[0] == -2

    def test_too_high_order(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_r = np.array([[2, 2, 0], [0, 1, 0]])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_status_3(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, 0, 0], [0, 0, 1]])
        X0 = np.array([10, 0, 0])
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(max_t=10, max_iter=100, chem_flag=True)
        assert sim.results.status_list[0] == 3

    def test_status_2(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        X0 = np.array([10, 0, 0])
        sim = Simulation(V_r, V_p, X0, k)
        sim.simulate(max_t=1, max_iter=100, chem_flag=True)
        assert sim.results.status_list[0] == 2

    def test_neg_k(self, setup_large):
        V_r, V_p, X0, k = setup_large
        k = np.array([1, 1, -1, 1, -1])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_Vp_Vr_shape(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, 1, 0]])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_kdet_Vr_shape(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        k = np.array([1, 1, 1])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_Vp_neg(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, -1, 0], [0, 0, 1]])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_Vr_neg(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        V_r = np.array([[-1, 0, 0], [0, 1, 0]])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_X0_neg(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        X0 = np.array([-10, 0, 0])
        with pytest.raises(ValueError):
            Simulation(V_r, V_p, X0, k)

    def test_reproduce(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(V_r, V_p, X0, k)
        sim2 = Simulation(V_r, V_p, X0, k)
        sim1.simulate()
        sim2.simulate()
        assert all(
            (i == j).all() for i, j in zip(sim1.results.t_list, sim2.results.t_list)
        )
        assert all(
            (i == j).all() for i, j in zip(sim1.results.x_list, sim2.results.x_list)
        )
        assert all(
            i == j for i, j in zip(sim1.results.status_list, sim2.results.status_list)
        )

    def test_reproduce_fail(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(V_r, V_p, X0, k)
        sim2 = Simulation(V_r, V_p, X0, k)
        sim1.simulate()
        sim2.simulate(seed=[1])
        assert not (
            all(
                (i == j).all() for i, j in zip(sim1.results.t_list, sim2.results.t_list)
            )
        )

    def test_incorrect_seed(self, setup_basic):
        V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(V_r, V_p, X0, k)
        with pytest.raises(ValueError):
            sim1.simulate(n_rep=2, seed=[1])
        with pytest.raises(ValueError):
            sim1.simulate(n_rep=2, seed=[1, 2, 3])


def test_bifurcation(setup_bifurcation):
    V_r, V_p, X0, k = setup_bifurcation
    count_excitation = 0
    n_runs = 100
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
