"""
    Tests for plotting function
"""


from cayenne.simulation import Simulation


def test_plotting(setup_basic):
    """ Test if plot generates without errors """
    species_name, rxn_names, V_r, V_p, X0, k = setup_basic
    n_runs = 10
    sim1 = Simulation(species_name, rxn_names, V_r, V_p, X0, k)
    sim1.simulate(
        algorithm="direct", max_t=150, max_iter=1000, chem_flag=True, n_rep=n_runs
    )
    sim1.plot(species_names=["A"])
    sim1.plot(species_names=["A", "B"], new_names=["Species A", "Secies B"])
    sim1.plot()
