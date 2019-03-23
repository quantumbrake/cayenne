pyssa : Python package for stochastic simulations
=================================================

|Build Status| |Updates| |Documentation Status| |pypi| |License| |Code
style: black|

Introduction
------------

``pyssa`` is a Python package for stochastic simulations. It offers a
simple api to define models, perform stochastic simulations on them and
visualize the results in a convenient manner.

Install
-------

Install with ``pip``:

.. code:: bash

   $ pip install pyssa

Documentation
-------------

-  General: https://pyssa.readthedocs.io.

Usage
-----

.. code:: python

   from pyssa.simulation import Simulation

    V_r = np.array([[1, 0], [0, 1], [0, 0]])  # Reactant matrix
    V_p = np.array([[0, 0], [1, 0], [0, 1]])  # Product matrix
    X0 = np.array([100, 0, 0])  # Initial state
    k = np.array([1.0, 1.0])  # Rate constants
    sim = Simulation(V_r, V_p, X0, k)  # Declare the simulation object
    # Run the simulation
    sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10)

You can change the algorithm used to perform the simulation by changing
the ``algorithm`` parameter

.. code:: python

   sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10, algorithm="tau_adaptive")

You can run the simulations on multiple cores by specifying the
``n_procs`` parameter

.. code:: python

   sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10, n_procs=4)

Plotting
~~~~~~~~

.. code:: python

   sim.plot()

.. image:: images/plot_basic.png

Plot of species A, B and C

Accessing the results
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   results = sim.results

Benchmarks
----------

We chose ``numba`` after extensive testing and benchmarking against
``python`` and ``cython`` implementations.

+---+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
| N | Min | Max | Mea | Std | Med | IQR | Out | OPS | Rou | Ite |
| a |     |     | n   | Dev | ian |     | lie |     | nds | rat |
| m |     |     |     |     |     |     | rs  |     |     | ion |
| e |     |     |     |     |     |     |     |     |     | s   |
| ( |     |     |     |     |     |     |     |     |     |     |
| t |     |     |     |     |     |     |     |     |     |     |
| i |     |     |     |     |     |     |     |     |     |     |
| m |     |     |     |     |     |     |     |     |     |     |
| e |     |     |     |     |     |     |     |     |     |     |
| i |     |     |     |     |     |     |     |     |     |     |
| n |     |     |     |     |     |     |     |     |     |     |
| m |     |     |     |     |     |     |     |     |     |     |
| s |     |     |     |     |     |     |     |     |     |     |
| ) |     |     |     |     |     |     |     |     |     |     |
+===+=====+=====+=====+=====+=====+=====+=====+=====+=====+=====+
| t | 314 | 342 | 322 | 11. | 318 | 9.1 | 1;1 | 3.0 | 5   | 1   |
| e | .17 | .99 | .93 | 459 | .79 | 533 |     | 966 |     |     |
| s | 58  | 15  | 18  | 0   | 83  | (1. |     | (1. |     |     |
| t | (1. | (1. | (1. | (1. | (1. | 0)  |     | 0)  |     |     |
| _ | 0)  | 0)  | 0)  | 0)  | 0)  |     |     |     |     |     |
| n |     |     |     |     |     |     |     |     |     |     |
| u |     |     |     |     |     |     |     |     |     |     |
| m |     |     |     |     |     |     |     |     |     |     |
| b |     |     |     |     |     |     |     |     |     |     |
| a |     |     |     |     |     |     |     |     |     |     |
| _ |     |     |     |     |     |     |     |     |     |     |
| b |     |     |     |     |     |     |     |     |     |     |
| e |     |     |     |     |     |     |     |     |     |     |
| n |     |     |     |     |     |     |     |     |     |     |
| c |     |     |     |     |     |     |     |     |     |     |
| h |     |     |     |     |     |     |     |     |     |     |
| m |     |     |     |     |     |     |     |     |     |     |
| a |     |     |     |     |     |     |     |     |     |     |
| r |     |     |     |     |     |     |     |     |     |     |
| k |     |     |     |     |     |     |     |     |     |     |
+---+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
| t | 17, | 19, | 18, | 862 | 18, | 1,0 | 2;0 | 0.0 | 5   | 1   |
| e | 345 | 628 | 255 | .47 | 148 | 30. |     | 548 |     |     |
| s | .76 | .39 | .39 | 11  | .93 | 367 |     | (0. |     |     |
| t | 98  | 31  | 31  | (75 | 58  | 6   |     | 02) |     |     |
| _ | (55 | (57 | (56 | .27 | (56 | (11 |     |     |     |     |
| c | .21 | .23 | .53 | )   | .93 | 2.5 |     |     |     |     |
| y | )   | )   | )   |     | )   | 7)  |     |     |     |     |
| _ |     |     |     |     |     |     |     |     |     |     |
| b |     |     |     |     |     |     |     |     |     |     |
| e |     |     |     |     |     |     |     |     |     |     |
| n |     |     |     |     |     |     |     |     |     |     |
| c |     |     |     |     |     |     |     |     |     |     |
| h |     |     |     |     |     |     |     |     |     |     |
| m |     |     |     |     |     |     |     |     |     |     |
| a |     |     |     |     |     |     |     |     |     |     |
| r |     |     |     |     |     |     |     |     |     |     |
| k |     |     |     |     |     |     |     |     |     |     |
+---+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
| t | 27, | 28, | 27, | 387 | 27, | 347 | 2;0 | 0.0 | 5   | 1   |
| e | 366 | 417 | 782 | .27 | 728 | .38 |     | 360 |     |     |
| s | .36 | .83 | .24 | 58  | .42 | 91  |     | (0. |     |     |
| t | 81  | 33  | 82  | (33 | 24  | (37 |     | 01) |     |     |
| _ | (87 | (82 | (86 | .80 | (86 | .95 |     |     |     |     |
| p | .11 | .85 | .03 | )   | .98 | )   |     |     |     |     |
| y | )   | )   | )   |     | )   |     |     |     |     |     |
| _ |     |     |     |     |     |     |     |     |     |     |
| b |     |     |     |     |     |     |     |     |     |     |
| e |     |     |     |     |     |     |     |     |     |     |
| n |     |     |     |     |     |     |     |     |     |     |
| c |     |     |     |     |     |     |     |     |     |     |
| h |     |     |     |     |     |     |     |     |     |     |
| m |     |     |     |     |     |     |     |     |     |     |
| a |     |     |     |     |     |     |     |     |     |     |
| r |     |     |     |     |     |     |     |     |     |     |
| k |     |     |     |     |     |     |     |     |     |     |
+---+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+

License
-------

Copyright (c) 2018-2019, Dileep Kishore, Srikiran Chandrasekaran. Released
under: Apache Software License 2.0

Credits
-------

This package was created with
`Cookiecutter <https://github.com/audreyr/cookiecutter>`__ and the
`audreyr/cookiecutter-pypackage <https://github.com/audreyr/cookiecutter-pypackage>`__
project template.

.. |Build Status| image:: https://travis-ci.com/Heuro-labs/pyssa.svg?token=qCMKydrUTvcJ87J6czex&branch=master
   :target: https://travis-ci.com/Heuro-labs/pyssa
.. |Updates| image:: https://pyup.io/repos/github/Heuro-labs/pyssa/shield.svg
   :target: https://pyup.io/repos/github/Heuro-labs/pyssa/
.. |Documentation Status| image:: https://readthedocs.org/projects/pyssa/badge/?version=latest
   :target: https://pyssa.readthedocs.io/en/latest/?badge=latest
.. |pypi| image:: https://img.shields.io/pypi/v/pyssa.svg
   :target: https://pypi.python.org/pypi/pyssa
.. |License| image:: https://img.shields.io/badge/license-Apache%202-blue.svg
.. |Code style: black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
