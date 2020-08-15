"""
    Contains various utility functions
"""

cimport cython
cimport numpy as np
import numpy as np
# import random


Na = 6.023e23  # Avogadro's constant
HIGH = 1e20
TINY = 1e-20


# @cython.returns((int, int))
@cython.boundscheck(False)
@cython.wraparound(False)
cdef roulette_selection(double [:] prop_list, long long [:] Xt):
    """
        Perform roulette selection on the list of propensities.

        Return the index of the selected reaction (``choice``) by performing
        Roulette selection on the given list of reaction propensities.

        Parameters
        ----------
        prop_list: array_like
            A 1D array of the propensities of the reactions.

        Xt: array_like
            A 1D array of the current simulation state.

        Returns
        -------
        choice: int
            Index of the chosen reaction.
        status: int
            Status of the simulation as described in ``direct``.
    """
    cdef:
        int choice= 0, status=0, counter=0, counter_max=prop_list.shape[0], Xt_counter_max = Xt.shape[0]
        double prop0 = 0.0, Xtsum = 0.0, TINY = 1e-20
    cdef (int, int) result
    for counter in range(counter_max):
        prop0 += prop_list[counter]
    if prop0 < TINY:
        for counter in range(Xt_counter_max):
            Xtsum += Xt[counter]
        counter = 0
        if Xtsum < TINY:
            status = 3
            result = -1, status
            return result
        else:
            status = -2
            result = -1, status
            return result

    cdef double prop_sum = prop_list[0]/prop0
    # r1 = random.random()  # Roll the wheel
    r1 = np.random.random()  # Roll the wheel
    for counter in range(1, counter_max + 1):
        if r1 < prop_sum:
            choice = counter - 1
            break
        else:
            prop_sum += prop_list[counter]/prop0
    result = choice, 0
    return result


def get_kstoc(react_stoic, k_det, volume, chem_flag):
    """
        Compute k_stoc from k_det.

        Return a vector of the stochastic rate constants (k_stoc) determined
        from the deterministic rate constants (k_det) [#]_

        Parameters
        ----------
        react_stoic: (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are columns and species are rows.
        k_det: (nr,) ndarray
            A 1D array representing the deterministic rate constants of the
            system.
        volume: float
            The volume of the reactor vessel which is important for second
            and higher order reactions
        chem_flag: bool
            If True, divide by Na (Avogadro's constant) while calculating
            stochastic rate constants. Defaults to ``False``.

        Returns
        -------
        k_stoc: (nr,) ndarray
            A 1D array representing the stochastic rate constants of the
            system.

        References
        ----------
        .. [#] Gillespie, D.T., 1976.
            A general method for numerically simulating the stochastic time evolution
            of coupled chemical reactions. J. Comput. Phys. 22, 403–434.
            doi:10.1016/0021-9991(76)90041-3.
    """
    cdef:
        int nr=react_stoic.shape[1], ns=react_stoic.shape[0], max_order
        double factor = 1.0
    orders = np.zeros((nr,), dtype=np.int64) # Order of rxn = number of reactants
    for ind1 in range(nr):
        for ind2 in range(ns):
            orders[ind1] += react_stoic[ind2, ind1]
    cdef double[:] k_stoc = k_det.copy()
    if chem_flag:
        factor = Na
    for ind in range(nr):
        max_order = max(react_stoic[:, ind])
        # If highest order is 3
        if max_order == 3:
            k_stoc[ind] = k_det[ind] * 6 / pow(factor * volume, 2)
        elif max_order == 2:  # Highest order is 2
            k_stoc[ind] = k_det[ind] * 2 / pow(factor * volume, orders[ind] - 1)
        else:
            k_stoc[ind] = k_det[ind] / pow(factor * volume, orders[ind] - 1)
    return k_stoc

def py_roulette_selection(prop_list, Xt):
    """
        Perform roulette selection on the list of propensities.

        Return the index of the selected reaction (``choice``) by performing
        Roulette selection on the given list of reaction propensities.

        Parameters
        ----------
        prop_list: array_like
            A 1D array of the propensities of the reactions.

        Xt: array_like
            A 1D array of the current simulation state.

        Returns
        -------
        choice: int
            Index of the chosen reaction.
        status: int
            Status of the simulation as described in ``direct``.
    """
    return roulette_selection(prop_list, Xt)
