# cython: profile=True

cimport cython
import numpy as np

Na = 6.023e23  # Avogadro's constant
HIGH = 1e20
TINY = 1e-20


def sumfunc(a, b):
    return a + b + b

@cython.boundscheck(False)
@cython.wraparound(False)
def roulette_selection(double[:] prop_list, int[:] Xt):
    """Perform roulette selection on the list of propensities.

    Return the index of the selected reaction (`choice`) by performing
    Roulette selection on the given list of reaction propensities.

    Parameters
    ----------
    prop_list : array_like
        A 1D array of the propensities of the reactions.

    Xt : array_like
        A 1D array of the current simulation state.

    Returns
    -------
    choice : int
        Index of the chosen reaction.
    status : int
        Status of the simulation as described in `direct`.
    """
    cdef int choice = 0
    cdef int status = 0
    cdef int counter = 0
    cdef int counter_max = prop_list.shape[0]
    cdef int Xt_counter_max = Xt.shape[0]
    cdef double prop0 = 0.0
    cdef double Xtsum = 0.0
    for counter in range(counter_max):
        prop0 += prop_list[counter]
    if prop0 < TINY:
        for counter in range(Xt_counter_max):
            Xtsum += Xt[counter]
        counter = 0
        if Xtsum < TINY:
            status = 3
            return -1, status
        else:
            status = -2
            return -1, status

    cdef double prop_sum = prop_list[0]/prop0
    r1 = np.random.rand()  # Roll the wheel
    for counter in range(1, counter_max + 1):
        if r1 < prop_sum:
            choice = counter - 1
            break
        else:
            prop_sum += prop_list[counter]/prop0
    return choice, 0


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
    cdef int nr = react_stoic.shape[1]
    orders = np.sum(react_stoic, axis=0)  # Order of rxn = number of reactants
    cdef float factor = 1.0
    k_stoc = k_det.copy()
    if chem_flag:
        factor = Na
    for ind in range(nr):
        # If highest order is 3
        if react_stoic[:, ind].max() == 3:
            k_stoc[ind] = k_det[ind] * 6 / np.power(factor * volume, 2)
        elif react_stoic[:, ind].max() == 2:  # Highest order is 2
            k_stoc[ind] = k_det[ind] * 2 / np.power(factor * volume, orders[ind] - 1)
        else:
            k_stoc[ind] = k_det[ind] / np.power(factor * volume, orders[ind] - 1)
    return k_stoc
