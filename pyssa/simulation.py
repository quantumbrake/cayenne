"""
    Naive implementation of the Gillespie algorithm (direct method) in Numba
"""

import numpy as np
from numba import njit
from typing import Tuple
from .direct_naive import direct_naive


class Simulation:
    """
        A main class for running simulations

        Parameters
        ---------
        react_stoic : (nr, ns) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are rows and species are columns.
        prod_stoic : (nr, ns) ndarray
            A 2D array of the stoichiometric coefficients of the products.
            Reactions are rows and species are columns.
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
        ---------
        results : Results
            The results instance

        Raises
        ------
        ValueError
            If supplied with order > 3.

        References
        ----------
        1. Gillespie, D.T., 1976. A general method for numerically
        simulating the stochastic time evolution of coupled chemical
        reactions. J. Comput. Phys. 22, 403–434.
        doi:10.1016/0021-9991(76)90041-3.
        2. Cao, Y., Gillespie, D.T., Petzold, L.R., 2006.
        Efficient step size selection for the tau-leaping simulation
        method. J. Chem. Phys. 124, 044109. doi:10.1063/1.2159468
        3. Gupta, A., 2013. Parameter estimation in deterministic
        and stochastic models of biological systems. University of
        Wisconsin-Madison.

        Examples
        --------
        >>> V_r = np.array([[1,0,0],[0,1,0]])
        >>> V_p = np.array([[0,1,0],[0,0,1]])
        >>> X0 = np.array([10,0,0])
        >>> k = np.array([1,1])
        >>> [_, _, status] = direct_naive(V_r, V_p, X0, k, max_t = 1, max_iter = 100)
    """

    results = None

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
        self._nr = self._react_stoic.shape[0]
        self._ns = self._react_stoic.shape[1]
        self._volume = volume
        self._orders = np.sum(
            self._react_stoic, 1
        )  # Order of rxn = number of reactants
        self._check_consistency()

    def _check_consistency(self):
        if (self._nr != self._prod_stoic.shape[0]) or (
            self._ns != self._prod_stoic.shape[1]
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

    def simulate(
        self,
        max_t: float = 10.0,
        max_iter: int = 1000,
        volume: float = 1.0,
        seed: int = 0,
        n_rep: int = 1,
        algorithm: str = "direct_naive",
        **kwargs,
    ):
        """
        Run the simulation

        Parameters
        ----------
        max_t : float, optional
            The end time of the simulation. The default is `max_t`=1.0 units.
        max_iter : int, optional
            The maximum number of iterations of the simulation loop. The
            default is 100 iterations.

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
        tlist = [0.0]
        xlist = []
        status_list = []
        if algorithm == "direct_naive":
            for _ in range(n_rep):
                t, X, status = direct_naive(
                    self._react_stoic,
                    self._prod_stoic,
                    self._init_state,
                    self._k_det,
                    max_t,
                    max_iter,
                    volume,
                    seed,
                    self._chem_flag,
                )
                tlist.append(t)
                xlist.append(X)
                status_list.append(status)
        else:
            raise ValueError("Requested algorithm not supported")

