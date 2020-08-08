.. figure:: https://raw.githubusercontent.com/Heuro-labs/cayenne/master/docs/images/logo.png
   :alt: Logo for cayenne

   Logo for cayenne

cayenne : Python package for stochastic simulations
===================================================

|Travis Build Status| |Azure Build Status| |codecov| |Updates|
|Documentation Status| |pypi| |License| |Code style: black|

Introduction
------------

``cayenne`` is a Python package for stochastic simulations. It offers a
simple API to define models, perform stochastic simulations with them
and visualize the results in a convenient manner.

Currently under active development in the ``develop`` branch.

Install
-------

Install with ``pip``:

.. code:: bash

   $ pip install cayenne

Documentation
-------------

-  General: https://cayenne.readthedocs.io.
-  Benchmark repository, comparing ``cayenne`` with other stochastic
   simulation packages: https://github.com/Heuro-labs/cayenne-benchmarks

Usage
-----

A short summary follows, but a more detailed tutorial can be found
`here <https://cayenne.readthedocs.io/en/latest/tutorial.html>`__. You
can define a model as a Python string (or a text file, see
`docs <https://cayenne.readthedocs.io>`__). The format of this string is
loosely based on the excellent
`antimony <https://tellurium.readthedocs.io/en/latest/antimony.html#introduction-basics>`__
library, which is used behind the scenes by ``cayenne``.

.. code:: python

   from cayenne.simulation import Simulation
   model_str = """
           const compartment comp1;
           comp1 = 1.0; # volume of compartment

           r1: A => B; k1;
           r2: B => C; k2;

           k1 = 0.11;
           k2 = 0.1;
           chem_flag = false;

           A = 100;
           B = 0;
           C = 0;
       """
   sim = Simulation.load_model(model_str, "ModelString")
   # Run the simulation
   sim.simulate(max_t=40, max_iter=1000, n_rep=10)
   sim.plot()

.. figure:: https://raw.githubusercontent.com/Heuro-labs/cayenne/master/docs/images/plot_basic.png
   :alt: Plot of species A, B and C

   Plot of species A, B and C

Change simulation algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can change the algorithm used to perform the simulation by changing
the ``algorithm`` parameter (one of ``"direct"``, ``"tau_leaping"`` or
``"tau_adaptive"``)

.. code:: python

   sim.simulate(max_t=150, max_iter=1000, n_rep=10, algorithm="tau_leaping")

Our `benchmarks <https://github.com/Heuro-labs/cayenne-benchmarks>`__
are summarized `below <#benchmarks>`__, and show ``direct`` to be a good
starting point. ``tau_leaping`` offers greater speed but needs
specification and tuning of the ``tau`` hyperparameter. The
``tau_adaptive`` is less accurate and a work in progress.

Run simulations in parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can run the simulations on multiple cores by specifying the
``n_procs`` parameter

.. code:: python

   sim.simulate(max_t=150, max_iter=1000, n_rep=10, n_procs=4)

Accessing simulation results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can access all the results or the results for a specific list of
species

.. code:: python

   # Get all the results
   results = sim.results
   # Get results only for one or more species
   results.get_species(["A", "C"])

You can also access the final states of all the simulation runs by

.. code:: python

   # Get results at the simulation endpoints
   final_times, final_states = results.final

Additionally, you can access the state a particular time point of
interest :math:`t`. ``cayenne`` will interpolate the value from nearby
time points to give an accurate estimate.

.. code:: python

   # Get results at timepoint "t"
   t = 10.0
   states = results.get_state(t) # returns a list of numpy arrays

.. raw:: html

   <h2 id="benchmarks">

Benchmarks

.. raw:: html

   </h2>

+-----------------+-----------------+-----------------+-----------------+
|                 | direct          | tau_leaping     | tau_adaptive    |
+=================+=================+=================+=================+
| cayenne         | :he             | :he             | Less accurate   |
|                 | avy_check_mark: | avy_check_mark: | than            |
|                 | Most accurate   | Very fast but   | GillespieSSA’s  |
|                 | yet             | may need manual | version         |
|                 |                 | tuning          |                 |
+-----------------+-----------------+-----------------+-----------------+
| Tellurium       | :exclamation:   | N/A             | N/A             |
|                 | Inaccurate for  |                 |                 |
|                 | 2nd order       |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| GillespieSSA    | Very slow       | :exclamation:   | :exclamation:   |
|                 |                 | Inaccurate for  | Inaccurate for  |
|                 |                 | initial zero    | initial zero    |
|                 |                 | counts          | counts          |
+-----------------+-----------------+-----------------+-----------------+
| BioSimulator.jl | :exclamation:   | :exclamation:   | :exclamation:   |
|                 | Inaccurate      | Inaccurate for  | Inaccurate for  |
|                 | interpolation   | initial zero    | initial zero    |
|                 |                 | counts          | counts          |
+-----------------+-----------------+-----------------+-----------------+

License
-------

Copyright (c) 2018-2020, Dileep Kishore, Srikiran Chandrasekaran.
Released under: Apache Software License 2.0

Credits
-------

-  `Cython <https://cython.org/>`__
-  `antimony <https://tellurium.readthedocs.io/en/latest/antimony.html>`__
-  `pytest <https://docs.pytest.org>`__
-  `Cookiecutter <https://github.com/audreyr/cookiecutter>`__
-  `audreyr/cookiecutter-pypackage <https://github.com/audreyr/cookiecutter-pypackage>`__
-  `black <https://github.com/ambv/black>`__
-  Logo made with `logomakr <https://logomakr.com/>`__

.. |Travis Build Status| image:: https://travis-ci.com/Heuro-labs/cayenne.svg?branch=master
   :target: https://travis-ci.com/Heuro-labs/cayenne
.. |Azure Build Status| image:: https://dev.azure.com/srikiranc/cayenne/_apis/build/status/Heuro-labs.cayenne?branchName=master
   :target: https://dev.azure.com/srikiranc/cayenne/_build
.. |codecov| image:: https://codecov.io/gh/Heuro-labs/cayenne/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/Heuro-labs/cayenne
.. |Updates| image:: https://pyup.io/repos/github/Heuro-labs/cayenne/shield.svg
   :target: https://pyup.io/repos/github/Heuro-labs/cayenne/
.. |Documentation Status| image:: https://readthedocs.org/projects/cayenne/badge/?version=latest
   :target: https://cayenne.readthedocs.io/en/latest/?badge=latest
.. |pypi| image:: https://img.shields.io/pypi/v/cayenne.svg
   :target: https://pypi.python.org/pypi/cayenne
.. |License| image:: https://img.shields.io/badge/license-Apache%202-blue.svg
.. |Code style: black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
