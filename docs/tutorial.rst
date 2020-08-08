Tutorial
========

Model Building
--------------

Consider a simple system of chemical reactions given by:

.. math::

    A \xrightarrow[]{k_1} B\\
    B \xrightarrow[]{k_2} C\\

Suppose k\ :sub:`0.11` = 1, k\ :sub:`0.1` = 1 and there are initially 100 units of A. Then we have the following model string ::

    >>> model_str = """
            const compartment comp1;
            comp1 = 1.0; # volume of compartment

            r1: A => B; k1; # differs from antimony
            r2: B => C; k2; # differs from antimony

            k1 = 0.11;
            k2 = 0.1;
            chem_flag = false;

            A = 100;
            B = 0;
            C = 0;
        """

The format of the model string is based on a subset of the `antimony modeling language <https://tellurium.readthedocs.io/en/latest/antimony.html#introduction-basics>`_, but with one key difference. ``Antimony`` allows the user to specify custom rate *equations* for each reaction. ``cayenne`` automagically generates the rate equations behind the scenes, and user need only supply the rate *constants*. The format is discussed below:

Model format
^^^^^^^^^^^^

::

    const compartment comp1;

This defines the compartment in which the reactions happen.::

    comp1 = 1.0;

This defines the volume of the compartment in which reactions happen. For zero and first order reactions, this number does not matter. For second and third order reactions, varying compartment volume will affect the kinetic outcomes even when the rest of the model is not changed. A blank line after this separates these definitions from the reactions.::

    r1: A => B; k1; # differs from antimony
    r2: B => C; k2; # differs from antimony

Here ``r1`` and ``r2`` refer to the names of the reactions. This is followed by a colon and the reactants in that reaction. In ``r1`` there is only one reactant, ``A``. Additional reactants or stoichiometries can be written like ``A + 2B``. This is followed by a ``=>`` which separates reactants and products. Products are written in a fashion similar to the reactants. A semi-colon indicates the end of the products. This is followed by a symbol depicting the rate constant e.g. ``k1``, and the reaction ends with a second semi-colon. A blank line after this separates these reactions from rate-constant assignments.::

    k1 = 0.11;
    k2 = 0.1;

The rate constants are assigned one per line, with each line ending in a semi-colon. Every rate constant defined in the reactions must be assigned a numerical value at this stage, or ``cayenne`` will throw a `cayenne.model_io.RateConstantError`.::

    chem_flag = false;

An additional element that is included at this stage is the ``chem_flag`` boolean variable. This is discussed more in detail in the documentation of cayenne.Simulation class under the notes section. Briefly, if

1. the system under consideration is a chemical system and the supplied rate constants are in units of molarity or M or mol/L, ``chem_flag`` should be set to ``true``
2. the system under consideration is a biological system and the supplied rate constants are in units of copies/L or CFU/L, ``chem_flag`` should be set to ``false``

A blank line after this separates rate constants from initial values for the species.::

    A = 100;
    B = 0;
    C = 0;

The initial values for species are assigned one per line, with each line ending in a semi-colon. Every species defined in the reactions must be assigned an integer initial value at this stage, or ``cayenne`` will throw a `cayenne.model_io.InitialStateError`.

.. warning::

    ``Antimony`` has a set of reserved keywords that cannot be used as species, compartment or variable names, Eg. ``formula``, ``priority``, ``time``, etc. Refer to the `antimony documentation <https://tellurium.readthedocs.io/en/latest/antimony.html#language-reference>`_ for more information.

.. note::

    ``cayenne`` only accepts zero, first, second and third order reactions. We decided to not allow custom rate equations for stochastic simulations for two reasons:

    1. A custom rate equation, such as the Monod equation (see here_ for background) equation below, may violate the assumptions_ of stochastic simulations. These assumptions include a well stirred chamber with molecules in Brownian motion, among others.

    .. math::

        \mu = \frac{\mu_{max}S}{K_S + S}

    2. An equation resembling the Monod equation, the Michaelis-Menten_ equation, is grounded chemical kinetic theory. Yet the rate expression (see below) does not fall under 0-3 order reactions supported by ``cayenne``. However, the *elementary* reactions that make up the Michaelis-Menten kinetics are first and second order in nature. These *elementary* reactions can easily be modeled with ``cayenne``, but with the specification of additional constants (see `examples <examples.html>`_). A study shows that using the rate expression of Michaelis-Menten kinetics is valid under `some conditions <https://pubmed.ncbi.nlm.nih.gov/21261403/>`_.


    .. math::

        \frac{dP}{dt} = \frac{\mu_{max}S}{K_S + S}

.. _antimony:
.. _here: https://en.wikipedia.org/wiki/Monod_equation
.. _assumptions: https://en.wikipedia.org/wiki/Gillespie_algorithm
.. _Michaelis-Menten: https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics

.. note::
    The ``chem_flag`` is set to ``True`` since we are dealing with a chemical system. For defintion of ``chem_flag``, see the notes under the definition of the ``Simulation`` class.

These variables are passed to the ``Simulation`` class to create an object that represents the current system ::

    >>> from cayenne import Simulation

    >>> sim = Simulation.load_model(model_str, "ModelString")

.. autoclass:: cayenne.simulation.Simulation

Running Simulations
-------------------

Suppose we want to run 10 repetitions of the system for at most 1000 steps / 40 time units each, we can use the ``simulate`` method to do this. ::

    >>> sim.simulate(max_t=40, max_iter=1000, n_rep=10)

.. automethod:: cayenne.simulation.Simulation.simulate


Plotting
--------

To plot the results on the screen, we can simply plot all species concentrations at all the time-points using: ::

    >>> sim.plot()

.. image:: ../docs/images/plot_basic.png
    :scale: 70%
    :align: center
    :alt: Plot of A, B and C species over time.


A subset of the species can be plotted along with custom display names by supplying additional arguments to ``Simulation.plot`` as follows:

    >>> sim.plot(species_names = ["A", "C"], new_names = ["Starting material", "Product"])

.. image:: ../docs/images/plot_ac.png
    :scale: 70%
    :align: center
    :alt: Plot of A and B species over time.

By default, calling the plot object returns the ``matplotlib`` figure and axis objects. To display the plot, we just do::

    >>> sim.plot()
    >>> import matplotlib.pyplot as plt
    >>> plt.show()

Instead to save the figure directly to a file, we do::

    >>> sim.plot()
    >>> import matplotlib.pyplot as plt
    >>> plt.savefig("plot.png")

.. automethod:: cayenne.simulation.Simulation.plot

.. note::

    1. The ``sim.plot`` method needs to be run after running ``sim.simulate``
    2. More detailed plots can be created manually by accessing the ``sim.results`` object


Accessing the results
---------------------

.. currentmodule:: cayenne.results

The results of the simulation can be retrieved by accessing the ``Results`` object as ::

    >>> results = sim.results
    >>> results
    <Results species=('A', 'B', 'C') n_rep=10algorithm=direct sim_seeds=[8325804 1484405 2215104 5157699 8222403 7644169 5853461 6739698  374564
 2832983]>

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

Additionally, you can also access a particular species' trajectory through time across all simulations with the ``get_species`` function as follows: ::

    >>> results.get_species(["A", "C"])

This will return a list with a numpy array for each repetition. We use a list here instead of higher dimensional ndarray for the following reason: any two repetitions of a stochastic simulation need not return the same number of time steps.

.. autoclass:: cayenne.results.Results

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

    >>> sim.simulate(max_t=150, max_iter=1000, n_rep=10, algorithm="tau_leaping")
