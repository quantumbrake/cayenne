"""
    The main class for running stochastic simulation
"""

from functools import partial
import multiprocessing as mp
from typing import List, Optional
from warnings import warn

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from .direct import direct
from .tau_leaping import tau_leaping
from .tau_adaptive import tau_adaptive
from .results import Results


def wrapper(x, func):
    return func(*x)


class Simulation:
    """
        A main class for running simulations.

        Parameters
        ----------
        react_stoic : (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are columns and species are rows.
        prod_stoic : (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the products.
            Reactions are columns and species are rows.
        init_state : (ns,) ndarray
            A 1D array representing the initial state of the system.
        k_det : (nr,) ndarray
            A 1D array representing the deterministic rate constants of the
            system.
        volume : float, optional
            The volume of the reactor vessel which is important for second
            and higher order reactions. Defaults to 1 arbitrary units.
        chem_flag : bool, optional
            If True, divide by Na while calculating stochastic rate constants.
            Defaults to False.

        Attributes
        ----------
        results : Results
            The results instance

        Raises
        ------
        ValueError
            If supplied with order > 3.

        Examples
        --------
        >>> V_r = np.array([[1,0],[0,1],[0,0]])
        >>> V_p = np.array([[0,0],[1,0],[0,1]])
        >>> X0 = np.array([10,0,0])
        >>> k = np.array([1,1])
        >>> sim = Simulation(V_r, V_p, X0, k)
        >>> sim.simulate(max_t=10, max_iter=100, n_rep=n_runs)
    """

    _results: Optional[Results] = None

    def __init__(
        self,
        react_stoic: np.ndarray,
        prod_stoic: np.ndarray,
        init_state: np.ndarray,
        k_det: np.ndarray,
        chem_flag: bool = False,
        volume: float = 1.0,
    ) -> None:
        self._react_stoic = react_stoic
        self._prod_stoic = prod_stoic
        self._init_state = init_state
        self._k_det = k_det
        self._chem_flag = chem_flag
        self._ns = self._react_stoic.shape[0]
        self._nr = self._react_stoic.shape[1]
        self._volume = volume
        self._orders = np.sum(
            self._react_stoic, axis=0
        )  # Order of rxn = number of reactants
        self._check_consistency()

    def _check_consistency(self):
        if (self._ns != self._prod_stoic.shape[0]) or (
            self._nr != self._prod_stoic.shape[1]
        ):
            raise ValueError("react_stoic and prod_stoic should be the same shape.")
        if np.any(self._react_stoic < 0):
            raise ValueError("react_stoic cannot have negative elements.")
        if np.any(self._prod_stoic < 0):
            raise ValueError("V_p cannot have negative elements.")
        if np.any(self._init_state < 0):
            raise ValueError("Initial numbers in X0 can't be negative.")
        if np.any(self._k_det < 0):
            raise ValueError("Rate constant(s) can't be negative.")
        if self._k_det.shape[0] != self._nr:
            raise ValueError("Number of rate constants must equal number of reactions")
        if self._chem_flag not in (True, False):
            raise ValueError("chem_flag must be a boolean True or False.")
        if np.max(self._orders) > 3:
            raise ValueError("Order greater than 3 not suppported.")

    @property
    def results(self) -> Optional[Results]:
        """
            The ``Results`` instance of the simulation

            Returns
            -------
            Optional[Results]
        """
        if self._results is None:
            warn(
                "Run `Simulation.simulate` before requesting the results object",
                Warning,
            )
            return self._results
        else:
            return self._results

    def simulate(
        self,
        max_t: float = 10.0,
        max_iter: int = 1000,
        volume: float = 1.0,
        seed: Optional[List[int]] = None,
        n_rep: int = 1,
        n_procs: int = 1,
        algorithm: str = "direct",
        **kwargs,
    ):
        """
        Run the simulation

        Parameters
        ----------
        max_t : float, optional
            The end time of the simulation
            The default is 10.0
        max_iter : int, optional
            The maximum number of iterations of the simulation loop
            The default is 1000 iterations
        volume : float, optional
            The volume of the system
            The default value is 1.0
        seed : List[int], optional
            The list of seeds for the simulations
            The length of this list should be equal to `n_rep`
            The default value is None
        n_rep : int, optional
            The number of repetitions of the simulation required
            The default value is 1
        n_procs : int, optional
            The number of cpu cores to use for the simulation
            The default value is 1
        algorithm : str, optional
            The algorithm to be used to run the simulation
            The default value is "direct"

        Returns
        -------
        t : float
            End time of the simulation.
        Xt : ndarray
            System state at time `t` and initial.
        status : int
            Indicates the status of the simulation at exit.
            1 : Succesful completion, terminated when `max_iter` iterations reached.
            2 : Succesful completion, terminated when `max_t` crossed.
            3 : Succesful completion, terminated when all species went extinct.
            -1 : Failure, order greater than 3 detected.
            -2 : Failure, propensity zero without extinction.
        """
        tlist = []
        xlist = []
        status_list = []

        if seed is not None:
            if n_rep != len(seed):
                raise ValueError("Seed should be as long as n_rep")
        else:
            seed = [index for index in range(n_rep)]

        if not isinstance(max_iter, int):
            raise TypeError("max_iter should be of type int")

        algo_args = []
        if algorithm == "direct":
            for index in range(n_rep):
                algo_args.append(
                    (
                        self._react_stoic,
                        self._prod_stoic,
                        self._init_state,
                        self._k_det,
                        max_t,
                        max_iter,
                        volume,
                        seed[index],
                        self._chem_flag,
                    )
                )
                algo = direct
        elif algorithm == "tau_leaping":
            if "tau" in kwargs.keys():
                tau = kwargs["tau"]
            else:
                tau = 0.1
            for index in range(n_rep):
                algo_args.append(
                    (
                        self._react_stoic,
                        self._prod_stoic,
                        self._init_state,
                        self._k_det,
                        tau,
                        max_t,
                        volume,
                        seed[index],
                        self._chem_flag,
                    )
                )
                algo = tau_leaping
        elif algorithm == "tau_adaptive":
            if "epsilon" in kwargs.keys():
                epsilon = kwargs["epsilon"]
            else:
                epsilon = 0.03
            if "nc" in kwargs.keys():
                nc = kwargs["nc"]
            else:
                nc = 10
            for index in range(n_rep):
                algo_args.append(
                    (
                        np.int64(self._react_stoic),
                        np.int64(self._prod_stoic),
                        self._init_state,
                        self._k_det,
                        nc,
                        epsilon,
                        max_t,
                        max_iter,
                        volume,
                        seed[index],
                        self._chem_flag,
                    )
                )
                algo = tau_adaptive
        else:
            raise ValueError("Requested algorithm not supported")
        algo_func = partial(wrapper, func=algo)
        with mp.Pool(processes=n_procs) as pool:
            results = pool.map(algo_func, algo_args)
            for t, X, status in results:
                tlist.append(t)
                xlist.append(X)
                status_list.append(status)
            self._results = Results(tlist, xlist, status_list, algorithm, seed)

    def plot(self, plot_indices: list = None, disp: bool = True, names: list = None):
        """
        Plot the simulation

        Parameters
        ----------
        plot_indices : list, optional
            The indices of the species to be plotted.
            The default is `[i for i in range(self._ns)]` plots all species.
        disp : bool, optional
            If `True`, the plot is displayed.
            The default shows the plot.
        names : list, optional
            The names of the species to be plotted.
            The default is `xi` for species `i`.

        Returns
        -------
        fig : class 'matplotlib.figure.Figure'
            Figure object of the generated plot.
        ax : class 'matplotlib.axes._subplots.AxesSubplot
            Axis objected of the generated plot.
        """

        if self._results is None:
            raise ValueError("Simulate not run.")
        else:
            if plot_indices is None:
                plot_indices = [i for i in range(self._ns)]
            elif np.any(np.array(plot_indices) < 0):
                raise ValueError("Negative indexing not supported")

            n_indices = len(plot_indices)
            prop_cycle = plt.rcParams["axes.prop_cycle"]
            colors = prop_cycle.by_key()["color"]
            fig, ax = plt.subplots()
            res = self._results
            legend_handlers = [0] * n_indices
            generic_names = [""] * n_indices
            for index1 in range(n_indices):
                legend_handlers[index1] = mlines.Line2D([], [], color=colors[index1])
                generic_names[index1] = "x" + str(plot_indices[index1])
                for index2 in range(len(res.status_list)):
                    ax.step(
                        res.t_list[index2],
                        res.x_list[index2][:, plot_indices[index1]],
                        color=colors[index1],
                    )
            if names is None:
                names = generic_names
            fig.legend(legend_handlers, names)
            if disp:
                plt.show()
            return fig, ax
