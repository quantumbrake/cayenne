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

`pyssa` is a Python package for stochastic simulations. It offers a simple api to define models, perform stochastic simulations on them and visualize the results in a convenient manner.


## Install

Install with `pip`:

```bash
$ pip install pyssa
```


## Documentation

  - General: <https://pyssa.readthedocs.io>.


## Usage

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

You can change the algorithm used to perform the simulation by changing the `algorithm` parameter

```python
sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10, algorithm="tau_adaptive")
```

You can run the simulations on multiple cores by specifying the `n_procs` parameter

```python
sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10, n_procs=4)
```


### Plotting

```python
sim.plot()
```

![Plot of species A, B and C](https://raw.githubusercontent.com/Heuro-labs/pyssa/master/docs/images/plot_basic.png)

### Accessing the results

```python
results = sim.results
```

You can also access the final states of all the simulation runs by

```python
final_times, final_states = results.final
```

## Benchmarks

We chose `numba` after extensive testing and benchmarking against `python` and `cython` implementations.

Name (time in ms) | Min | Max | Mean | StdDev | Median | IQR | Outliers | OPS | Rounds | Iterations |
| --- |:---:| :---:|:---:| :---:|:---:| :---:|:---:| :---:|:---:| :---:|
test_numba_benchmark | 314.1758 (1.0) | 342.9915 (1.0) | 322.9318 (1.0) |  11.4590 (1.0) |   318.7983 (1.0) |  9.1533 (1.0) |     1;1 | 3.0966 (1.0) |     5 |     1 |
test_cy_benchmark |  17,345.7698 (55.21)  |  19,628.3931 (57.23)  |  18,255.3931 (56.53)  | 862.4711 (75.27) |  18,148.9358 (56.93)  |  1,030.3676 (112.57) |   2;0 | 0.0548 (0.02) |    5 |     1 |
test_py_benchmark |  27,366.3681 (87.11) |  28,417.8333 (82.85) |   27,782.2482 (86.03)  |  387.2758 (33.80)  |  27,728.4224 (86.98)  |  347.3891 (37.95) |   2;0 | 0.0360 (0.01) |    5 |     1 |

## License

Copyright (c) 2018-2019, Dileep Kishore, Srikiran Chandrasekaran. Released under: Apache Software License 2.0

## Credits

- [Numba](https://numba.pydata.org/)
- [pytest](https://docs.pytest.org)
- [Cookiecutter](https://github.com/audreyr/cookiecutter)
- [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage)
- [black](https://github.com/ambv/black)
