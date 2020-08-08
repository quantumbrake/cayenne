=======
History
=======

1.0.2 (2020-08-08)
------------------

Added
+++++
- Warnings about antimony keyword usage in tutorial, ``ModelIO`` class

Fixed
+++++
- ``setup.py`` setup requirements are automatically installed
- ``antimony`` version incompatibility issue

1.0.0 (2020-07-18)
------------------

Added
+++++
- ``Antimony`` support and ``ModelIO`` class, giving easier entry point to load models
- Support for custom species names in plotting and ``Results``
- Support for automatic ``cpu`` core detection
- New logo for ``cayenne``

Changed
+++++++
- Package name changed from ``pyssa`` to ``cayenne``
- Update docs for new API

0.9.1 (2020-05-23)
------------------
Fix interpolation bug in ``Results.get_state``

Fixed
+++++
- ``Results.get_state`` function now adds an epsilon to time

0.9.0 (2019-12-14)
------------------
Replace ``numba`` implementation with ``Cython`` implementation

Fixed
+++++
- Propensity calculation in all algorithms
- Considerable speed-up in algorithm runtimes
- Remove ``volume`` from ``Simulation.simulate`` parameters

Added
+++++
- ``direct`` algorithm in ``Cython``
- ``tau_leaping`` algorithm in ``Cython``
- ``tau_adaptive`` algorithm in ``Cython`` (experimental)
- ``Cython`` to Azure pipeline
- Accuracy tests from `sbml-test-suite <https://github.com/sbmlteam/sbml-test-suite>`_
- HOR property and tests for it
- Code coverage for ``Cython``
- Algorithms page to the documentation
- Examples page to the documentation

Changed
+++++++
- Remove ``numba`` algorithms
- Remove interpolation for direct algorithm
- ``sim.plot`` now plots ``post`` step curve
- Updated tutorial page of the documentation


0.8.2 (2019-04-20)
------------------

Fixed
+++++
- Initialize ``algorithms`` submodule with ``__init__.py``
- Update ``setup.py`` to allow submodule detection

0.8.0 (2019-04-13)
------------------

Added
+++++
- ``Results.get_states`` method - returns state at time ``t``
- Accuracy tests for all algorithms
- Additional consistency checks for ``X0`` and ``k_det``

Changed
+++++++
- Refactor algorithms into sub module ``algorithms``
- Refactor algorithm independent tests

Fixed
+++++
- Indexing issue in propensity calculation in ``direct`` algorithm
- Indexing issue in propensity calculation in ``tau_leaping`` algorithm
- Address edge case X->2X in ``tau_adaptive`` algorithm

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
