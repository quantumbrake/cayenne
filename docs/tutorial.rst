Tutorial
========

Model Building
--------------

Consider a simple system of chemical reactions given by:

.. math::

    A \xrightarrow[]{k_1} B\\
    B \xrightarrow[]{k_2} C\\

Suppose k\ :sub:`1` = 1, k\ :sub:`2` = 1 and there are initially 100 units of A. Then we have the following variable definitions ::

    >>> k1, k2 = 1.0, 1.0
    >>> A0, B0, C0 = 100, 0, 0

Then to build the model we have the following variable definitions::

    >>> import numpy as np

    # reaction stoichiometry
    >>> V_r = np.array([[1, 0], [0, 1], [0, 0]])
    # product stoichiometry
    >>> V_p = np.array([[0, 0], [1, 0], [0, 1]])
    # initial concentration
    >>> X0 = np.array([A0, B0, C0])
    # reaction rates
    >>> k = np.array([k1, k2])

We then pass these variables to the ``Simulation`` class to create an object that represents the current system ::

    >>> from pyssa import Simulation

    >>> sim = Simulation(V_r, V_p, X0, k)

.. autoclass:: pyssa.simulation.Simulation

Running Simulations
-------------------

Suppose we want to run 10 repetitions of the system for at most 1000 steps / 150 time units each, we can use the ``simulate`` method to do this. ::

    >>> from pyssa import Simulation

    >>> sim = Simulation(V_r, V_p, X0, k)
    >>> sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10)

.. automethod:: pyssa.simulation.Simulation.simulate

.. note::
    The ``chem_flag`` is set to ``True`` since we are dealing with a chemical system.

Plotting
--------

To plot the results on the screen, we can simply plot all species concentrations at all the time-points using: ::

    >>> sim.plot()

.. image:: ../docs/images/plot_basic.png
    :scale: 70%
    :align: center
    :alt: Plot of A, B and C species over time.

To plot only certain species concentrations we can use the ``plot_indices`` parameter. These correspond to the indices of the species assigned during the initialization of the ``Simulation`` instance. For example to plot only A and B, we use - ``[0,1]`` ::

    >>> sim.plot(plot_indices = [0, 1])

.. image:: ../docs/images/plot_ab.png
    :scale: 70%
    :align: center
    :alt: Plot of A and B species over time.

To not display the plot on the screen and retrieve the figure and axis objects, we set the ``disp`` parameter to False ::

    >>> fig, ax = sim.plot(disp = False)

.. automethod:: pyssa.simulation.Simulation.plot

.. note::

    1. The ``sim.plot`` method needs to be run after running ``sim.simulate``
    2. More detailed plots can be created manually by accessing the ``sim.results`` object


Accessing the results
---------------------

.. currentmodule:: pyssa.results

The results of the simulation can be retrieved by accessing the ``Results`` object as ::

    >>> results = sim.results
    >>> results
    <Results n_rep=10 algorithm=direct sim_seeds=[8325804 1484405 2215104 5157699 8222403 7644169 5853461 6739698  374564 2832983]>

The ``Results`` object provides abstractions for easy retrieval and iteration over the simulation results. For example you can iterate over every run of the simulation using ::

    >>> for x, t, status in results:
    ...     pass

You can access the results of the ``n`` th run by ::

    >>> nth_result = results[n]

You can also access the final states of all the simulation runs by ::

    >>> final_times, final_states = results.final

    # final times of each repetition
    >>> final_times
    array([6.23502469, 7.67449057, 6.15181435, 8.95810706, 7.12055223,
       7.06535004, 6.07045973, 7.67547689, 9.4218006 , 9.00615099])

    # final states of each repetition
    >>> final_states
    array([[  0,   0, 100],
       [  0,   0, 100],
       [  0,   0, 100],
       [  0,   0, 100],
       [  0,   0, 100],
       [  0,   0, 100],
       [  0,   0, 100],
       [  0,   0, 100],
       [  0,   0, 100],
       [  0,   0, 100]])

You can obtain the state of the system at a particular time using the ``get_state`` method. For example to get the state of the system at time ``t=5.0`` for each repetition: ::

    >>> results.get_state(5.0)
    [array([ 1,  4, 95]),
    array([ 1,  2, 97]),
    array([ 0,  2, 98]),
    array([ 3,  4, 93]),
    array([ 0,  3, 97]),
    array([ 0,  2, 98]),
    array([ 1,  1, 98]),
    array([ 0,  4, 96]),
    array([ 1,  6, 93]),
    array([ 1,  3, 96])]

.. autoclass:: pyssa.results.Results

.. autosummary::

    Results.__iter__
    Results.__len__
    Results.__contains__
    Results.__getitem__
    Results.final
    Results.get_state

Algorithms
----------

The ``Simulation`` class currently supports the following algorithms (see :ref:`algorithms`):

1. :ref:`Gillespie's direct method <direct>`
2. :ref:`Tau leaping method method <tau_leaping>`
3. :ref:`Adaptive tau leaping method (experimental) <tau_adaptive>`

You can change the algorithm used to perform a simulation using the ``algorithm`` argument ::

    >>> sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10, algorithm="tau_leaping")
