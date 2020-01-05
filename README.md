# pyssa : Python package for stochastic simulations

[![Build Status](https://travis-ci.com/Heuro-labs/pyssa.svg?token=qCMKydrUTvcJ87J6czex&branch=master)](https://travis-ci.com/Heuro-labs/pyssa)
[![Build Status](https://dev.azure.com/srikiranc/pyssa/_apis/build/status/Heuro-labs.pyssa?branchName=master)](https://dev.azure.com/srikiranc/pyssa/_build/latest?definitionId=1?branchName=master)
[![codecov](https://img.shields.io/codecov/c/github/Heuro-labs/pyssa.svg)](https://codecov.io/gh/Heuro-labs/pyssa)
[![Updates](https://pyup.io/repos/github/Heuro-labs/pyssa/shield.svg)](https://pyup.io/repos/github/Heuro-labs/pyssa/)
[![Documentation Status](https://readthedocs.org/projects/pyssa/badge/?version=latest)](https://pyssa.readthedocs.io/en/latest/?badge=latest)
[![pypi](https://img.shields.io/pypi/v/pyssa.svg)](https://pypi.python.org/pypi/pyssa)
![License](https://img.shields.io/badge/license-Apache%202-blue.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)



## Introduction

`pyssa` is a Python package for stochastic simulations. It offers a simple API to define models, perform stochastic simulations with them and visualize the results in a convenient manner.

Currently under active development in the `develop` branch.

## Install

Install with `pip`:

```bash
$ pip install pyssa
```


## Documentation

  - General: <https://pyssa.readthedocs.io>.


## Usage

A short summary follows, but a more detailed tutorial can be found at <https://pyssa.readthedocs.io/en/latest/tutorial.html>

```python
from pyssa.simulation import Simulation
V_r = np.array([[1, 0], [0, 1], [0, 0]])  # Reactant matrix
V_p = np.array([[0, 0], [1, 0], [0, 1]])  # Product matrix
X0 = np.array([100, 0, 0])  # Initial state
k = np.array([1.0, 1.0])  # Rate constants
sim = Simulation(V_r, V_p, X0, k)  # Declare the simulation object
# Run the simulation
sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10)
```

### Change simulation algorithm

You can change the algorithm used to perform the simulation by changing the `algorithm` parameter

```python
sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10, algorithm="tau_adaptive")
```

### Run simulations in parallel
You can run the simulations on multiple cores by specifying the `n_procs` parameter

```python
sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10, n_procs=4)
```

### Plot simulation results

```python
sim.plot()
```

![Plot of species A, B and C](https://raw.githubusercontent.com/Heuro-labs/pyssa/master/docs/images/plot_basic.png)

### Accessing simulation results

```python
results = sim.results
```

You can also access the final states of all the simulation runs by

```python
final_times, final_states = results.final
```

## License

Copyright (c) 2018-2020, Dileep Kishore, Srikiran Chandrasekaran. Released under: Apache Software License 2.0

## Credits

- [Cython](https://cython.org/)
- [pytest](https://docs.pytest.org)
- [Cookiecutter](https://github.com/audreyr/cookiecutter)
- [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage)
- [black](https://github.com/ambv/black)
