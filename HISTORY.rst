=======
History
=======

0.5.4 (2018-12-02)
------------------

Added
+++++
- badge to readme

0.5.3 (2018-12-02)
------------------

Added
+++++
- plot to pypi

Changed
+++++++
- fix bumpversion/black issue
- remove history from package long_description


0.5.0 (2018-12-01)
------------------

First public release!!

Added
+++++
- testpypi deployment
- pyup security checking
- readthedocs deployment
- Tutorials and documentation
- Plotting functionality through ``Simulation.plot``

Changed
+++++++
- ``Simulation.results`` is now a property
- Updated tests to support the new api changes

Chore
+++++
- Updated the README


0.4.0 (2018-11-23)
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
