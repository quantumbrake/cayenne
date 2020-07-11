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

A short summary follows, but a more detailed tutorial can be found [here](https://pyssa.readthedocs.io/en/latest/tutorial.html). You can define a model as a Python string (or a text file, see [docs](https://pyssa.readthedocs.io)). The format of this string is loosely based on the excellent [antimony](https://tellurium.readthedocs.io/en/latest/antimony.html#introduction-basics) library, which is used behind the scenes by `pyssa`.

```python
from pyssa.simulation import Simulation
import matplotlib.pyplot as plt
model_str = """
        const compartment comp1;
        comp1 = 7; # volume of compartment

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
```

![Plot of species A, B and C](https://raw.githubusercontent.com/Heuro-labs/pyssa/master/docs/images/plot_basic.png)


### Change simulation algorithm

You can change the algorithm used to perform the simulation by changing the `algorithm` parameter

```python
sim.simulate(max_t=150, max_iter=1000, n_rep=10, algorithm="tau_adaptive")
```

### Run simulations in parallel
You can run the simulations on multiple cores by specifying the `n_procs` parameter

```python
sim.simulate(max_t=150, max_iter=1000, n_rep=10, n_procs=4)
```

### Accessing simulation results

You can access all the results or the results for a specific list of species

```python
# Get all the results
results = sim.results
# Get results only for one or more species
results.get_species(["A", "C"])
```

You can also access the final states of all the simulation runs by

```python
# Get results at the simulation endpoints
final_times, final_states = results.final
```

Additionally, you can access the state a particular time point of interest $t$. `pyssa` will interpolate the value from nearby time points to give an accurate estimate.

```python
# Get results at timepoint "t"
t = 10.0
states = results.get_state(t) # returns a list of numpy arrays
```

## License

Copyright (c) 2018-2020, Dileep Kishore, Srikiran Chandrasekaran. Released under: Apache Software License 2.0

## Credits

- [Cython](https://cython.org/)
- [antimony](https://tellurium.readthedocs.io/en/latest/antimony.html)
- [pytest](https://docs.pytest.org)
- [Cookiecutter](https://github.com/audreyr/cookiecutter)
- [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage)
- [black](https://github.com/ambv/black)
