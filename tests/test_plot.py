"""
    Tests for plotting function
"""


from pyssa.simulation import Simulation
import pytest


def test_plotting(setup_basic):
    V_r, V_p, X0, k = setup_basic
    n_runs = 10
    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=n_runs)
    sim1.plot(plot_indices=[0, 1], disp=False)
    assert 1 == 1
    with pytest.raises(ValueError):
        sim1.plot(plot_indices=[-1], disp=False)
    sim1.plot(plot_indices=[0, 1], names=["Species A", "Species B"], disp=False)
