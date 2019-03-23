"""
    Module containing some utility functions
"""

from typing import Tuple

import numpy as np
from numba import njit

Na = 6.023e23  # Avogadro's constant
HIGH = 1e20
TINY = 1e-20


@njit(nogil=True, cache=False)
def get_kstoc(
    react_stoic: np.ndarray, k_det: np.ndarray, volume: float, chem_flag: bool
) -> np.ndarray:
    """Compute k_stoc from k_det.

    Return a vector of the stochastic rate constants (k_stoc) determined
    from the deterministic rate constants (k_det).

    Parameters
    ----------
    react_stoic : (ns, nr) ndarray
        A 2D array of the stoichiometric coefficients of the reactants.
        Reactions are columns and species are rows.
    k_det : (nr,) ndarray
        A 1D array representing the deterministic rate constants of the
        system.
    volume : float
        The volume of the reactor vessel which is important for second
        and higher order reactions
    chem_flag : bool
        If True, divide by Na while calculating stochastic rate constants.

    Returns
    -------
    k_stoc : (nr,) ndarray
        A 1D array representing the stochastic rate constants of the
        system.

    References
    ----------
    1. Gillespie, D.T., 1976. A general method for numerically
    simulating the stochastic time evolution of coupled chemical
    reactions. J. Comput. Phys. 22, 403–434.
    doi:10.1016/0021-9991(76)90041-3.
    """
    nr = react_stoic.shape[1]
    orders = np.sum(react_stoic, axis=0)  # Order of rxn = number of reactants
    k_stoc = k_det.copy()
    if chem_flag:
        factor = Na
    else:
        factor = 1.0
    for ind in range(nr):
        # If highest order is 3
        if react_stoic[:, ind].max() == 3:
            k_stoc[ind] = k_det[ind] * 6 / np.power(factor * volume, 2)
        elif react_stoic[:, ind].max() == 2:  # Highest order is 2
            k_stoc[ind] = k_det[ind] * 2 / np.power(factor * volume, orders[ind] - 1)
        else:
            k_stoc[ind] = k_det[ind] / np.power(factor * volume, orders[ind] - 1)
    return k_stoc


@njit(nogil=True, cache=False)
def roulette_selection(prop_list: np.ndarray, Xt: np.ndarray) -> Tuple[int, int]:
    """Perform roulette selection on the list of propensities.

    Return the index of the selected reaction (`choice`) by performing
    Roulette selection on the given list of reaction propensities.

    Parameters
    ----------
    prop : array_like
        A 1D array of the propensities of the reactions.

    Returns
    -------
    choice : int
        Index of the chosen reaction.
    status : int
        Status of the simulation as described in `direct`.
    """
    prop0 = np.sum(prop_list)  # Sum of propensities
    # choice = 0
    if prop0 < TINY:
        if np.sum(Xt) < TINY:
            status = 3
            return -1, status
        else:
            status = -2
            return -1, status
    prop_norm = prop_list / prop0  # Normalize propensities to be < 1
    # Concatenate 0 to list of probabilities
    probs = np.cumsum(prop_norm)
    r1 = np.random.rand()  # Roll the wheel
    # Identify where it lands and update that reaction VERIFY MAY BE WRONG
    choice = np.searchsorted(probs, r1)
    return choice, 0
