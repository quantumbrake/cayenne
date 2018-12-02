# Stochastic Simulation Algorithms in Python

[![Build Status](https://travis-ci.com/Heuro-labs/pyssa.svg?token=qCMKydrUTvcJ87J6czex&branch=master)](https://travis-ci.com/Heuro-labs/pyssa)
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

sim = Simulation(V_r, V_p, X0, k)
sim.simulate(max_t=150, max_iter=1000, chem_flag=True, n_rep=10)
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

## Benchmarks

We chose `numba` after extensive testing and benchmarking against `python` and `cython` implementations.

Name (time in ms) | Min | Max | Mean | StdDev | Median | IQR | Outliers | OPS | Rounds | Iterations |
| --- |:---:| :---:|:---:| :---:|:---:| :---:|:---:| :---:|:---:| :---:|
test_numba_benchmark | 314.1758 (1.0) | 342.9915 (1.0) | 322.9318 (1.0) |  11.4590 (1.0) |   318.7983 (1.0) |  9.1533 (1.0) |     1;1 | 3.0966 (1.0) |     5 |     1 |
test_cy_benchmark |  17,345.7698 (55.21)  |  19,628.3931 (57.23)  |  18,255.3931 (56.53)  | 862.4711 (75.27) |  18,148.9358 (56.93)  |  1,030.3676 (112.57) |   2;0 | 0.0548 (0.02) |    5 |     1 |
test_py_benchmark |  27,366.3681 (87.11) |  28,417.8333 (82.85) |   27,782.2482 (86.03)  |  387.2758 (33.80)  |  27,728.4224 (86.98)  |  347.3891 (37.95) |   2;0 | 0.0360 (0.01) |    5 |     1 |

## License

Copyright (c) 2018, Dileep Kishore, Srikiran Chandrasekaran. Released under: Apache Software License 2.0

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage) project template.
