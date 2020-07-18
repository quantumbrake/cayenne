"""
    Implementation of the
    `tau leaping algorithm <https://en.wikipedia.org/wiki/Tau-leaping>`_.
    This is an approximate method that needs to be tuned to the system at hand
    (by modifying the time step or the ``tau`` parameter).
    A default ``tau=0.1`` is assumed by
    ``cayenne``. This algorithm is approximate and faster than the Direct
    algorithm, but it must be used with caution. Smaller time steps make the
    simulation more accurate, but increase the code run time. Larger time steps
    makes the simulations less accurate and speeds up code run time.
"""

cimport cython
cimport numpy as np
import numpy as np
import random
from ..utils import get_kstoc, TINY


@cython.boundscheck(False)
@cython.wraparound(False)
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
):
    """
        Runs the Tau Leaping Simulation Algorithm.
        Exits if negative population encountered.

        Parameters
        ----------
        react_stoic: (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are columns and species are rows.
        prod_stoic: (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the products.
            Reactions are columns and species are rows.
        init_state: (ns,) ndarray
            A 1D array representing the initial state of the system.
        k_det: (nr,) ndarray
            A 1D array representing the deterministic rate constants of the
            system.
        tau: float
            The constant time step used to tau leaping.
        max_t: float
            The maximum simulation time to run the simulation for.
        volume: float
            The volume of the reactor vessel which is important for second
            and higher order reactions. Defaults to 1 arbitrary units.
        seed: int
            The seed for the numpy random generator used for the current run
            of the algorithm.
        chem_flag: bool
            If True, divide by Na (Avogadro's constant) while calculating
            stochastic rate constants. Defaults to ``False``.

        Returns
        -------
        t: ndarray
            Numpy array of the times.
        x: ndarray
            Numpy array of the states of the system at times in in `t`.
        status: int
            Indicates the status of the simulation at exit.

            1: Succesful completion, terminated when `max_iter` iterations reached.

            2: Succesful completion, terminated when `max_t` crossed.

            3: Succesful completion, terminated when all species went extinct.

            -1: Failure, order greater than 3 detected.

            -2: Failure, propensity zero without extinction.

            -3: Negative species count encountered.
    """
    cdef:
        int ite = 1, max_iter=0, ind1=0, ind2=0, xtsum=0  # Iteration counter
        double t_curr = 0.0, prop_sum=0.0
        Py_ssize_t ns=react_stoic.shape[0], nr=react_stoic.shape[1]
    v = prod_stoic - react_stoic  # ns x nr
    xt = init_state.copy().astype(np.int64)  # Number of species at time t_curr
    max_iter = np.int(max_t / tau) + 1
    x = np.zeros((max_iter, ns), dtype=np.int64)
    t = np.zeros((max_iter))
    x[0, :] = init_state.copy()
    n_events = np.zeros((nr,), dtype=np.int64)
    np.random.seed(seed)  # Set the seed
    # Determine kstoc from kdet and the highest order or reactions
    prop = get_kstoc(react_stoic, k_det, volume, chem_flag)  # Vector of propensities
    kstoc = prop.copy()  # Stochastic rate constants

    cdef double [:] kstoc_view = kstoc
    cdef double [:] prop_view = prop
    cdef long long [:] xt_view = xt
    cdef long long [:] n_events_view = n_events
    cdef long [:, :] v_view = v
    cdef long [:, :] react_stoic_view = react_stoic
    cdef long long [:, :] x_view = x

    for ind in range(nr):
        prop_sum += prop_view[ind]
    if prop_sum < TINY:
        for ind in range(ns):
            xtsum += xt_view[ind]
        if xtsum > TINY:
            status = -2
            return t[:ite], x[:ite, :], status

    while ite < max_iter:
        # Calculate propensities
        prop_sum = 0
        for ind1 in range(nr):
            for ind2 in range(ns):
                # prop = kstoc * product of (number raised to order)
                if react_stoic_view[ind2, ind1]:
                    if react_stoic_view[ind2, ind1] == 1:
                        prop_view[ind1] *= x_view[ite - 1, ind2]
                    elif react_stoic_view[ind2, ind1] == 2:
                        prop_view[ind1] *= x_view[ite - 1, ind2] * (x_view[ite - 1, ind2] - 1) / 2
                    elif react_stoic_view[ind2, ind1] == 3:
                        prop_view[ind1] *= x_view[ite - 1, ind2] * (x_view[ite - 1, ind2] - 1) * (x_view[ite - 1, ind2] - 2) / 6
            prop_sum += prop_view[ind1]
            n_events_view[ind1] = np.random.poisson(prop_view[ind1] * tau)  # 1 x nr
        if prop_sum == 0:
            status = 3
            return t[:ite], x[:ite,], status
        for ind1 in range(ns):
            xt_view[ind1] = x_view[ite-1, ind1]
            for ind2 in range(nr):
                xt_view[ind1] += n_events_view[ind2] * v_view[ind1, ind2]
        for ind1 in range(ns):
            if xt_view[ind1] < 0:
                return t[:ite], x[:ite, :], -3
            else:
                x_view[ite, ind1] = xt_view[ind1]
        for ind1 in range(nr):
            prop_view[ind1] = kstoc_view[ind1]
        t_curr += tau
        t[ite] = t_curr
        ite += 1
    status = 2
    return t[:ite], x[:ite, :], status
