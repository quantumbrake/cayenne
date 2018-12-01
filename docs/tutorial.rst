==============
Model Building
==============


Consider a simple system of chemical reactions given by:

.. math::

    A \xrightarrow[]{k_1} B\\
    B \xrightarrow[]{k_2} C \\

Suppose k\ :sub:`1` = 1, k\ :sub:`2` 1 and there are initiall 100 units of A. Then we have the following variable definitions ::

    k1, k2 = 1.0, 1.0
    A0, B0, C0 = 100, 0, 0

Then to build the model we have the following variable definitions::

    import numpy as np
    V_r = np.array([[1, 0, 0], [0, 1, 0]])
    V_p = np.array([[0, 1, 0], [0, 0, 1]])
    X0 = np.array([A0, B0, C0])
    k = np.array([k1, k2])


===================
Running Simulations
===================

Suppose we want to run 10 runs of the system for earlier of 1000 time steps / 150 time units each, we have ::

    from pyssa.simulation import Simulation

    sim1 = Simulation(V_r, V_p, X0, k)
    sim1.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10)

Note that the ``chem_flag`` is set to ``True`` since we are dealing with a chemical system.

========
Plotting
========

To plot the results on the screen, we simply have ::

    sim1.plot()

To plot only A and B, we use the species indices (``[0,1]``) ::

    sim1.plot(plot_indices = [0, 1])

To not display the plot on the screen and retrieve the figure and axis objects, we have ::

    fig, ax = sim1.plot(disp = False)


=====================
Accessing the results
=====================

The results of the simulation can be retrieved by accessing the ``Results`` object as ::

    res = sim.results()
