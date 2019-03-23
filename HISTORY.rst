=======
History
=======

0.7.1 (2019-03-23)
------------------

Changed
+++++++
- Refactor ``tau_adaptive``
- Rename ``direct_naive`` to ``direct``

Fixed
+++++
- SSA part of ``tau_adaptive``
- Bug in linux compatibility of ``tau_adaptive``

0.7.0 (2019-02-02)
------------------

Added
+++++
- Support for the ``tau_adaptive`` algorithm
- Support for multiprocessing

Fixed
+++++
- Transpose stoichiometric matrix
- Update references in docstrings

Changed
+++++++
- Use ``TINY`` and ``HIGH`` for status estimation
- Use ``np.int64`` and ``np.float64`` explicitly

Chore
+++++
- Update dependencies
- Add azure pipelines for testing on Windows

0.6.0 (2018-12-16)
------------------

Added
+++++
- Updated ``direct_naive`` docstring
- Support for the ``tau_leaping`` algorithm
- Species name support for plotting

Fixed
+++++
- Check for sum propensities uses threshold instead of equality
- Add check for type of ``max_iter``

Changed
+++++++
- Update ``roulette_selection`` to use `np.searchsorted`
- Minor changes to ``numpy`` style usage

Chore
+++++
- Add ``codecov``
- Travis pypi autodepolyment
- Parameterize tests with algorithm name
- Add details about ``tau_leaping`` to docs and README


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
