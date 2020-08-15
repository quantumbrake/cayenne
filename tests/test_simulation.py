"""
    Tests for direct function
"""

import numpy as np
import pytest

from cayenne.simulation import Simulation


@pytest.mark.parametrize("algorithm", ["direct", "tau_leaping", "tau_adaptive"])
@pytest.mark.usefixtures("setup_basic", "setup_large")
class TestSanitizeAlg:
    """
        Sanity checks on Simulation class where simulations are attempted.
    """

    def test_null(self, algorithm, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        k = np.array([0.0, 0.0])
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm)
        assert sim.results.status_list[0] == -2

    def test_status_3(self, algorithm, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, 0], [0, 0], [0, 1]])
        X0 = np.array([10, 0, 0], dtype=np.int64)
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, max_t=10, max_iter=100, chem_flag=True)
        assert sim.results.status_list[0] == 3

    def test_status_2(self, algorithm, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        X0 = np.array([10, 0, 0], dtype=np.int64)
        sim = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim.simulate(algorithm=algorithm, max_t=1, max_iter=100, chem_flag=True)
        assert sim.results.status_list[0] == 2

    def test_reproduce(self, algorithm, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim2 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
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
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim2 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim1.simulate(algorithm=algorithm)
        sim2.simulate(algorithm=algorithm, seed=1)
        for i, j in zip(sim1.results.t_list, sim2.results.t_list):
            if i.shape[0] == j.shape[0]:
                assert not (
                    all(
                        (i == j).all()
                        for i, j in zip(sim1.results.t_list, sim2.results.t_list)
                    )
                )

    def test_incorrect_seed(self, algorithm, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        with pytest.raises(TypeError):
            sim1.simulate(algorithm=algorithm, n_rep=2, seed=[1])
        with pytest.raises(TypeError):
            sim1.simulate(algorithm=algorithm, n_rep=2, seed=[1, 2, 3])

    def test_maxiter_type(self, algorithm, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        with pytest.raises(TypeError):
            sim1.simulate(algorithm=algorithm, max_iter=100.0)

    def test_monotonic(self, algorithm, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        sim1 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim1.simulate(algorithm=algorithm)
        results = sim1.results
        for t_array in results.t_list:
            assert (np.diff(t_array) > 0).all()

    def test_long(self, algorithm, setup_long):
        species_names, rxn_names, V_r, V_p, X0, k = setup_long
        sim1 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        sim1.simulate(
            algorithm=algorithm, max_t=1e5, max_iter=int(1e8), chem_flag=False
        )
        X_output = np.array([0, 0, X0.sum()])
        assert (sim1.results.final[-1] == X_output).all()


@pytest.mark.usefixtures("setup_basic", "setup_large")
class TestSanitize:
    """
        Sanity checks on Simulation class where instance creation fails.
    """

    def test_too_high_order(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        V_r = np.array([[2, 0], [2, 1], [0, 0]])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_neg_k(self, setup_large):
        species_names, rxn_names, V_r, V_p, X0, k = setup_large
        k = np.array([1, 1, -1, 1, -1])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_k_shape(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        k = np.array([[1], [0], [0]])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_Vp_Vr_shape(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0], [1], [0]])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_kdet_Vr_shape(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        k = np.array([1, 1, 1])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_X0_Vr_shape(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        X0 = np.array([100, 0, 0, 0])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        X0 = np.array([100, 0])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_X0_2d(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        X0 = np.array([[100, 0, 0]])
        print(X0.shape)
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)
        X0 = np.array([[100], [0], [0]])
        print(X0.shape)
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_Vp_neg(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        V_p = np.array([[0, 0], [-1, 0], [0, 1]])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_Vr_neg(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        V_r = np.array([[-1, 0], [0, 1], [0, 0]])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)

    def test_X0_neg(self, setup_basic):
        species_names, rxn_names, V_r, V_p, X0, k = setup_basic
        X0 = np.array([-10, 0, 0])
        with pytest.raises(ValueError):
            Simulation(species_names, rxn_names, V_r, V_p, X0, k)


# def test_bifurcation(setup_bifurcation):
#     V_r, V_p, X0, k = setup_bifurcation
#     count_excitation = 0
#     n_runs = 10
#     deviation_tolerance = 0.05
#     sim1 = Simulation(species_names, rxn_names, V_r, V_p, X0, k)
#     sim1.simulate(
#         algorithm="tau_adaptive", max_t=150, max_iter=1000, chem_flag=True, n_rep=n_runs
#     )
#     for ind in range(n_runs):
#         x = sim1.results.x_list[ind]
#         xt = x[-1, :]
#         assert np.all(xt - np.array([0, 11, 0, 0]) == 0) or np.all(
#             xt - np.array([0, 0, 1, 10]) == 0
#         )
#         if np.all(xt - np.array([0, 11, 0, 0]) == 0):
#             count_excitation += 1
#     assert np.abs(count_excitation / n_runs - 0.5) < deviation_tolerance


class TestHOR:
    """
        Test the HOR method for various combinations of reactants.
    """

    @staticmethod
    def create_sim_inst(react_stoic):
        prod_stoic = np.random.randint(0, 5, react_stoic.shape)
        k = np.random.rand(react_stoic.shape[1])
        X = np.random.randint(0, 5, react_stoic.shape[0], dtype=np.int64)
        species_names = ["X"] * react_stoic.shape[0]
        rxn_names = ["R"] * react_stoic.shape[1]
        sim = Simulation(species_names, rxn_names, react_stoic, prod_stoic, X, k)
        return sim

    def test_hor1(self):
        # A ->
        react_stoic = np.array([[1, 0]])
        sim = self.create_sim_inst(react_stoic)
        assert sim.HOR == 1

    def test_hor2(self):
        # A ->
        # A ->
        # A ->
        react_stoic = np.array([[1, 1, 1], [0, 0, 0]])
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [1, 0])

    def test_hor3(self):
        # A + B ->
        # A ->
        react_stoic = np.array([[1, 1], [1, 0]])
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [2, 2])

    def test_hor4(self):
        # # A + B ->
        # # C ->
        react_stoic = np.array([[1, 0], [1, 0], [0, 1]])
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [2, 2, 1])

    def test_hor5(self):
        # A + A ->
        # A ->
        # B ->
        react_stoic = np.array([[2, 1, 0], [0, 0, 1]])
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [-2, 1])

    def test_hor6(self):
        # # A + A + B ->
        # # A + B ->
        # # B ->
        react_stoic = np.array([[2, 1, 0], [1, 1, 1]])
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [-32, 3])

    def test_hor7(self):
        # A + A + B ->
        # A + B + B ->
        react_stoic = np.array([[2, 1], [1, 2]])
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [-32, -32])

    def test_hor8(self):
        # # A + A + A ->
        # # A + A + B ->
        # # A + B ->
        react_stoic = np.array([[3, 2, 1], [0, 1, 1]])
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [-3, 3])

    def test_hor9(self):
        # A + A + A ->
        # B + B + B ->
        # A + A + B ->
        # A + B + B ->
        # A + B + C ->
        react_stoic = np.array([[3, 0, 2, 1, 1], [0, 3, 1, 2, 1], [0, 0, 0, 0, 1]])
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [-3, -3, 3])

    def test_hor10(self):
        # A + B ->
        # A + A + A ->
        # B + B ->
        # A + C ->
        # C + D + E ->
        react_stoic = np.array(
            [
                [1, 3, 0, 1, 0],
                [1, 0, 2, 0, 0],
                [0, 0, 0, 1, 1],
                [0, 0, 0, 0, 1],
                [0, 0, 0, 0, 1],
            ]
        )
        sim = self.create_sim_inst(react_stoic)
        assert np.all(sim.HOR == [-3, -2, 3, 3, 3])


class TestLoadModel:
    def test_correct_model_str(self, setup_00001_correct):
        sim = Simulation.load_model(setup_00001_correct, "ModelString")
        sim.simulate()
        assert sim.results.status_list == [2]
        sim.simulate(max_t=20.0, max_iter=5)
        assert sim.results.status_list == [1]

    def test_correct_model_file(self):
        sim = Simulation.load_model("tests/models/00001.txt", "ModelFile")
        sim.simulate()
        assert sim.results.status_list == [2]
