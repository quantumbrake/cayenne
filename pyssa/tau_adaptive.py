"""
    Implementation of the tau leaping algorithm in Numba
"""

from typing import Tuple
from numba import njit
import numpy as np
from .utils import get_kstoc

HIGH = 1e20


def HOR(react_stoic: np.ndarray):
    """
    Determine the HOR vector. HOR(i) is the highest order of reaction
    in which species S_i appears as a reactant.
        Parameters
        ---------
        react_stoic : (n_s, n_r) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are rows and species are columns.
    """
    n_s = react_stoic.shape[0]
    HOR = np.zeros([n_s])
    orders = np.sum(react_stoic, axis=0)
    for ind in range(n_s):
        this_orders = orders[np.where(react_stoic[ind, :] > 0)]
        if len(this_orders) == 0:
            HOR[ind] = 0
            continue
        HOR[ind] = np.max(this_orders)
        if HOR[ind] == 1:
            continue
        order_2_indices = np.where(orders == 2)
        if order_2_indices[0].size > 0:
            if np.max(react_stoic[ind, np.where(orders == 2)]) == 2 and HOR[ind] == 2:
                HOR[ind] = -2  # g_i should be (2 + 1/(x_i-1))
        if np.where(orders == 3):
            if (
                HOR[ind] == 3
                and np.max(react_stoic[ind, np.where(this_orders == 3)]) == 2
            ):
                HOR[ind] = -32  # g_i should be (3/2 * (2 + 1/(x_i-1)))
            elif (
                HOR[ind] == 3
                and np.max(react_stoic[ind, np.where(this_orders == 3)]) == 3
            ):
                HOR[ind] = -3  # g_i should be(3 + 1/(x_i-1) + 2/(x_i-2))
    return HOR


# @njit(nogil=True, cache=False)
def tau_adaptive(
    react_stoic: np.ndarray,
    prod_stoic: np.ndarray,
    init_state: np.ndarray,
    k_det: np.ndarray,
    nc: int,
    eps: float,
    max_t: float,
    volume: float,
    seed: int,
    chem_flag: bool,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
        Parameters
        ---------
        react_stoic : (n_r, n_s) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are rows and species are columns.
        prod_stoic : (n_r, n_s) ndarray
            A 2D array of the stoichiometric coefficients of the products.
            Reactions are rows and species are columns.
        init_state : (n_s,) ndarray
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
    epsilon = 0.03
    ite = 1  # Iteration counter
    t_curr = 0.0  # Time in seconds
    n_r = react_stoic.shape[0]
    n_s = react_stoic.shape[1]
    prod_stoic = np.transpose(prod_stoic)
    react_stoic = np.transpose(react_stoic)
    v = prod_stoic - react_stoic  # n_s x n_r
    xt = init_state.copy()  # Number of species at time t_curr
    max_iter = 10
    x = np.zeros((max_iter, n_s))
    t = np.zeros((max_iter))
    x[0, :] = init_state.copy()
    n_events = np.zeros((n_r,), dtype=np.int32)
    np.random.seed(seed)  # Set the seed
    # Determine kstoc from kdet and the highest order or reactions
    prop = np.copy(
        get_kstoc(react_stoic, k_det, volume, chem_flag)
    )  # Vector of propensities
    kstoc = prop.copy()  # Stochastic rate constants

    # Determine the HOR vector. HOR(i) is the highest order of reaction
    # in which species S_i appears as a reactant.
    HOR = np.zeros([n_s])
    orders = np.sum(react_stoic, axis=0)
    for ind in range(n_s):
        this_orders = orders[np.where(react_stoic[ind, :] > 0)]
        HOR[ind] = np.max(this_orders)
        if react_stoic[ind, np.argmax(this_orders)] == 2 and HOR[ind] == 2:
            HOR[ind] = -2  # g_i should be (2 + 1/(x_i-1))
        elif react_stoic[ind, np.argmax(this_orders)] == 2 and HOR[ind] == 3:
            HOR[ind] = -32  # g_i should be (3/2 * (2 + 1/(x_i-1)))
        elif react_stoic[ind, np.argmax(this_orders)] == 3:
            HOR[ind] = -3

    if np.sum(prop) < 1e-30:
        if np.sum(xt) > 1e-30:
            status = -2
            return t[:ite], x[:ite, :], status

    M = n_r
    N = n_s
    L = np.zeros(M)
    vis = np.zeros(M)
    react_species = np.where(np.sum(react_stoic, axis=1) > 0)[0]
    n_react_species = react_species.shape[0]
    mup = np.zeros(N)
    sigp = np.zeros(N)

    while ite < max_iter:
        # 1. Determine critical reactions

        # Calculate the propensities
        for ind1 in range(n_r):
            for ind2 in range(n_s):
                # prop = kstoc * product of (number raised to order)
                prop[ind1] *= np.power(xt[ind2], react_stoic[ind2, ind1])

        for ind in range(M):
            vis = react_stoic[:, ind]
            L[ind] = np.nanmin(x[ite - 1, :] / vis)
        crit = (L < nc) * (prop > 0)
        not_crit = True ^ crit
        # 2. Generate candidate taup
        if np.sum(not_crit) == 0:
            taup = HIGH
        else:
            # Compute mu from eqn 32a and sig from eqn 32b
            for ind, species_index in enumerate(react_species):
                temp = v[species_index, not_crit] * prop[not_crit]
                mup[ind] = np.sum(temp)
                sigp[ind] = np.sum(v[species_index, not_crit] * temp)
                if HOR[species_index] > 0:
                    g = HOR[species_index]
                elif HOR[species_index] == -2:
                    g = 1 + 2 / (xt[species_index] - 1)
                mup_list[ind] = max(epsilon * xt[species_index], 1)

            taup = 2
        # 3. For small taup, do SSA

        # 4. Generate second candidate taupp

        # 5. Leap

        # 6. Handle negatives

        # Calculate the event rates
        for ind1 in range(n_r):
            for ind2 in range(n_s):
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


if __name__ == "__main__":
    V_r = np.array([[1, 0, 0], [0, 1, 0]])
    V_p = np.array([[0, 1, 0], [0, 0, 1]])
    X0 = np.array([100, 0, 0])
    k = np.array([1.0, 1.0])
    tau_adaptive(
        V_r, V_p, X0, k, nc=10, eps=0.03, max_t=10, volume=1, seed=0, chem_flag=False
    )

