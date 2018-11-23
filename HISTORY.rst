=======
History
=======

0.3.0 (2018-11-23)
------------------

Added
+++++
- ``Simulation`` class - main class for running simulations
- ``Results`` class - for storing and acessing simulation results
- ``Simulation.simulate`` function that returns an instance of the ``Results`` class

Changed
+++++++
- Refactor ``get_kstoc`` and ``roulette_selection`` into ``utils.py``
- Refactor ``direct_naive`` into ``direct_naive.py``
- Delete ``pyssa.py`` and replace with ``Simulation`` class

Chore
+++++
- Add license and code-style badges
- Use ``black`` for code-formatting


0.2.0 (2018-11-10)
------------------

Added
+++++

- Naive implementation of the Gillepsie algorithm in ``numba``
- Tests - sanity checks, bifurcation and long running simulation
- CI on ``travis``


0.1.0 (2018-08-08)
------------------

* First commit
