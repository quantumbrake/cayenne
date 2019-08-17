"""
    Accuracy tests of the algorithms
"""
import numpy as np
import pytest

from pyssa.simulation import Simulation


def calculate_zy(sim: Simulation, time_list: list, mu_list: list, std_list: list):
    z_list = []
    y_list = []
    n_rep = len(sim.results.t_list)
    for ind1, t in enumerate(time_list[1:]):
        results = sim.results.get_state(t)
        mu_obs = np.mean(results)
        std_obs = np.std(results)
        z_list.append(
            np.sqrt(n_rep) * (mu_obs - mu_list[ind1 + 1]) / std_list[ind1 + 1]
        )
        y_list.append(
            np.sqrt(n_rep / 2) * ((std_obs ** 2) / (std_list[ind1 + 1] ** 2) - 1)
        )
    return np.array(z_list), np.array(y_list)


def calculate_zy_2sp(sim: Simulation, time_list: list, mu_list: list, std_list: list):
    z_list = []
    y_list = []
    n_rep = len(sim.results.t_list)
    for ind1, t in enumerate(time_list[1:]):
        results = sim.results.get_state(t)
        mu_obs = np.mean(results, axis=0)
        std_obs = np.std(results, axis=0)
        z_list.append(
            np.sqrt(n_rep) * (mu_obs - mu_list[ind1 + 1]) / std_list[ind1 + 1]
        )
        y_list.append(
            np.sqrt(n_rep / 2) * ((std_obs ** 2) / (std_list[ind1 + 1] ** 2) - 1)
        )
    return np.array(z_list), np.array(y_list)


@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00001(setup_00001, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00001
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00003(setup_00003, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00003
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    t_final, x_final = sim1.results.final
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert (np.array(sim1.results.status_list) != 1).all()
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00004(setup_00004, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00004
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    t_final, x_final = sim1.results.final
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert (np.array(sim1.results.status_list) != 1).all()
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00005(setup_00005, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00005
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    t_final, x_final = sim1.results.final
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert (np.array(sim1.results.status_list) != 1).all()
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00011(setup_00011, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00011
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00020(setup_00020, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00020
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00021(setup_00021, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00021
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00022(setup_00022, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00022
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00023(setup_00023, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00023
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00030(setup_00030, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00030
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy_2sp(sim1, time_list, mu_list, std_list)
    assert (np.array(sim1.results.status_list) != 1).all()
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00031(setup_00031, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00031
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy_2sp(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00037(setup_00037, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00037
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


@pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
def test_00038(setup_00038, algorithm):
    V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00038
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(
        algorithm=algorithm,
        max_t=max_t,
        max_iter=max_iter,
        chem_flag=False,
        n_rep=n_rep,
    )
    Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
    assert np.all(sim1.results.final[0] > time_list[-1])
    assert (-3 < Z).all()
    assert (Z < 3).all()
    assert (-5 < Y).all()
    assert (Y < 5).all()


# @pytest.mark.filterwarnings("ignore::UserWarning")
# @pytest.mark.parametrize("algorithm", ["direct_cython", "tau_leaping_cython"])
# def test_00039(setup_00039, algorithm):
#     V_r, V_p, X0, k, time_list, mu_list, std_list, max_t, max_iter, n_rep = setup_00039
#     sim1 = Simulation(V_r, V_p, X0, k)
#     sim1.simulate(
#         algorithm=algorithm,
#         max_t=max_t,
#         max_iter=max_iter,
#         chem_flag=False,
#         n_rep=n_rep,
#     )
#     print(sim1.results.status_list)
#     sim1.plot()
#     Z, Y = calculate_zy(sim1, time_list, mu_list, std_list)
#     assert np.all(sim1.results.final[0] > time_list[-1])
#     assert (-3 < Z).all()
#     assert (Z < 3).all()
#     assert (-5 < Y).all()
#     assert (Y < 5).all()
