"""
    Implementation of the tau leaping algorithm in Numba
"""

from typing import Tuple

from numba import njit
import numpy as np

from .utils import get_kstoc


@njit(nogil=True, cache=False)
def tau_leaping(
    react_stoic: np.ndarray,
    prod_stoic: np.ndarray,
    init_state: np.ndarray,
    k_det: np.ndarray,
    tau: float,
    max_t: float,
    volume: float,
    seed: int,
    chem_flag: bool,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
        Parameters
        ---------
        react_stoic : (n_r, ns) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are rows and species are columns.
        prod_stoic : (n_r, ns) ndarray
            A 2D array of the stoichiometric coefficients of the products.
            Reactions are rows and species are columns.
        init_state : (ns,) ndarray
            A 1D array representing the initial state of the system.
        k_det : (n_r,) ndarray
            A 1D array representing the deterministic rate constants of the
            system.
        tau : float
            The constant time step used to tau leaping.
        max_t : float
            The maximum simulation time to run the simulation for.
        volume : float
            The volume of the reactor vessel which is important for second
            and higher order reactions. Defaults to 1 arbitrary units.
        seed : int
            The seed for the numpy random generator used for the current run
            of the algorithm.
        chem_flag : bool
            If True, divide by Na while calculating stochastic rate constants.
            Defaults to False.

        Returns
        -------
        t : ndarray
            Numpy array of the times.
        x : ndarray
            Numpy array of the states of the system at times in in `t`.
        status : int
            Indicates the status of the simulation at exit.
            1 : Succesful completion, terminated when `max_iter` iterations reached.
            2 : Succesful completion, terminated when `max_t` crossed.
            3 : Succesful completion, terminated when all species went extinct.
            -1 : Failure, order greater than 3 detected.
            -2 : Failure, propensity zero without extinction.
            -3 : Negative species count encountered
    """

    ite = 1  # Iteration counter
    t_curr = 0.0  # Time in seconds
    n_r = react_stoic.shape[0]
    ns = react_stoic.shape[1]
    v = prod_stoic - react_stoic  # n_r x ns
    xt = init_state.copy()  # Number of species at time t_curr
    max_iter = np.int(max_t / tau) + 1
    x = np.zeros((max_iter, ns))
    t = np.zeros((max_iter))
    x[0, :] = init_state.copy()
    n_events = np.zeros((n_r,), dtype=np.int32)
    np.random.seed(seed)  # Set the seed
    # Determine kstoc from kdet and the highest order or reactions
    prop = np.copy(
        get_kstoc(react_stoic, k_det, volume, chem_flag)
    )  # Vector of propensities
    kstoc = prop.copy()  # Stochastic rate constants

    if np.sum(prop) < 1e-30:
        if np.sum(xt) > 1e-30:
            status = -2
            return t[:ite], x[:ite, :], status

    while ite < max_iter:
        # Calculate the event rates
        for ind1 in range(n_r):
            for ind2 in range(ns):
                # prop = kstoc * product of (number raised to order)
                prop[ind1] *= np.power(xt[ind2], react_stoic[ind1, ind2])
            n_events[ind1] = np.random.poisson(prop[ind1] * tau)  # 1 x n_r
        if np.all(prop == 0):
            status = 3
            return t[:ite], x[:ite,], status
        for ind1 in range(n_r):
            xt += n_events[ind1] * v[ind1, :]
        if np.any(xt < 0):
            return t[:ite], x[:ite, :], -3
        prop = np.copy(kstoc)
        x[ite, :] = xt
        t[ite] = t_curr
        t_curr += tau
        ite += 1
    status = 2
    return t[:ite], x[:ite, :], status
