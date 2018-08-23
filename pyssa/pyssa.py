"""
    Naive implementation of the Gillespie algorithm (Direct method)
"""

from typing import Tuple

import numpy as np


def direct_naive(
        V_r: np.ndarray,
        V_p: np.ndarray,
        X0: np.ndarray,
        k_det, max_t: float = 1.0,
        max_iter: int = 100,
        volume: float = 1.0,
        seed: int = 0,
) -> Tuple[float, np.ndarray, int]:
    """Naive implementation of the Direct method.

    A naive implementation of the Direct method of the SSA algorithm.

    Parameters
    ----------
    V_r : (nr, ns) ndarray
        A 2D array of the stoichiometric coefficients of the reactants.
        Reactions are rows and species are columns.
    V_p : (nr, ns) ndarray
        A 2D array of the stoichiometric coefficients of the products.
        Reactions are rows and species are columns.
    X0 : (ns,) ndarray
        A 1D array representing the initial state of the system.
    k_det : (nr,) ndarray
        A 1D array representing the deterministic rate constants of the
        system.
    max_t : float, optional
        The end time of the simulation. The default is `max_t`=1.0 units.
    max_iter : int, optional
        The maximum number of iterations of the simulation loop. The
        default is 100 iterations.
    volume : float, optional
        The volume of the reactor vessel which is important for second
        and higher order reactions. Defaults to 1 arbitrary units.

    Returns
    -------
    t : float
        End time of the simulation.
    Xt : ndarray
        System state at time `t` and initial.
    status : int
        Indicates the status of the simulation at exit.
            1 : Succesful completion, terminated when `max_iter`
            iterations reached.
            2 : Succesful completion, terminated when `max_t` croosed.
            3 : Succesful completion, terminated when all species
            went extinct.
            -1 : Failure, order greater than 3 detected.
            -2 : Failure, propensity zero without extinction.

    Raises
    ------
    ValueError
        If supplied with order > 3.

    References
    ----------

    .. [1] Gillespie, D.T., 1976. A general method for numerically
    simulating the stochastic time evolution of coupled chemical
    reactions. J. Comput. Phys. 22, 403–434.
    doi:10.1016/0021-9991(76)90041-3.
    .. [2] Cao, Y., Gillespie, D.T., Petzold, L.R., 2006.
    Efficient step size selection for the tau-leaping simulation
    method. J. Chem. Phys. 124, 044109. doi:10.1063/1.2159468
    .. [3] Gupta, A., 2013. Parameter estimation in deterministic
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

    ite = 1  # Iteration counter
    t = 0  # Time in seconds
    nr = V_r.shape[0]  # Number of reactions
    ns = V_r.shape[1]  # Number of species

    if (nr != V_p.shape[0]) or (ns != V_p.shape[1]):
        raise ValueError('V_r and V_p should be the same shape.')

    if (nr != k_det.shape[0]):
        raise ValueError('Number of elements in k_det must equal\
         number of rows in V_r.')

    if np.any(V_r < 0):
        raise ValueError('V_r cannot have negative elements.')

    if np.any(V_p < 0):
        raise ValueError('V_p cannot have negative elements.')

    if np.any(X0 < 0):
        raise ValueError('Initial numbers in X0 can\'t be negative.')

    if np.any(k_det < 0):
        neg_indices = np.where(k_det < 0)[0]
        raise ValueError('Rate constant(s) at position(s) '+str(neg_indices)+' are negative.')

    V = V_p - V_r  # nr x ns
    Xt = np.copy(X0)  # Number of species at time t
    Xtemp = np.zeros(nr)  # Temporary X for updating
    kstoc = np.zeros(nr)  # Stochastic rate constants
    orders = np.sum(V_r, 1)  # Order of rxn = number of reactants
    status = 0
    np.random.seed(seed = seed)  # Set the seed

    if np.max(orders) > 3:
        raise ValueError('Order greater than 3 detected.')

    if np.max(orders) > 1:
        raise RuntimeWarning('Order greater than 1, using volume = ', volume)

    # Determine kstoc from kdet and the highest order or reactions
    kstoc = get_kstoc(k_det, V_r, volume)
    prop = np.copy(kstoc) # Vector of propensities

    while ite < max_iter:
        # Calculate propensities
        for ind1 in range(nr):
            for ind2 in range(ns):
                # prop = kstoc * product of (number raised to order)
                prop[ind1] *= np.power(Xt[ind2], V_r[ind1, ind2])
        # Roulette wheel
        [choice, status] = roulette_selection(prop,Xt)
        if status == 0:
            Xtemp = Xt + V[choice, :]
        else:
            return t, Xt, status
        print(Xt, Xtemp)
        print('-' * 80)
        # If negative species produced, reject step
        if np.min(Xtemp) < 0:
            continue
        # Update Xt and t
        else:
            Xt = Xtemp
            r2 = np.random.rand()
            t += 1 / np.sum(prop) * np.log(1 / r2)
            if t > max_t:
                status = 2
                print("Reached maximum time (t = )",t)
                return t, Xt, status
        prop = np.copy(kstoc)
        ite += 1
    status = 1
    return t, Xt, status

def roulette_selection(prop_list, Xt):
    r"""Perform roulette selection on the list of propensities.

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
        Status of the simulation as described in `direct_naive`.

    Examples
    --------
    >>>

    """
    prop = np.copy(prop_list)
    prop0 = np.sum(prop)  # Sum of propensities
    if prop0 == 0:
        if np.sum(Xt) == 0:
            status = 3
            return [-1, status]
        else:
            status = -2
            return [-1, status]
    prop = prop / prop0  # Normalize propensities to be < 1
    # Concatenate 0 to list of probabilities
    probs = [0]+list(np.cumsum(prop))
    r1 = np.random.rand() # Roll the wheel
    # Identify where it lands and update that reaction
    for ind1 in range(len(probs)):
        if r1 <= probs[ind1]:
            choice = ind1 - 1
            break
    return [choice, 0]

def get_kstoc(k_det, V_r, volume = 1.0):
    r"""Compute k_stoc from k_det.

    Return a vector of the stochastic rate constants (k_stoc) determined
    from the deterministic rate constants (k_det).

    Parameters
    ----------
    k_det : (nr,) ndarray
        A 1D array representing the deterministic rate constants of the
        system.
    V_r : (nr, ns) ndarray
        A 2D array of the stoichiometric coefficients of the reactants.
        Reactions are rows and species are columns.
    volume : float, optional
        The volume of the reactor vessel which is important for second
        and higher order reactions. Defaults to 1 arbitrary units.

    Returns
    -------
    k_stoc : (nr,) ndarray
        A 1D array representing the stochastic rate constants of the
        system.

    References
    ----------

    .. [1] Gillespie, D.T., 1976. A general method for numerically
    simulating the stochastic time evolution of coupled chemical
    reactions. J. Comput. Phys. 22, 403–434.
    doi:10.1016/0021-9991(76)90041-3.

    Examples
    --------
    >>>

    """

    Na = 6.023e23  # Avogadro's constant
    nr = V_r.shape[0]  # Number of reactions
    orders = np.sum(V_r, 1)  # Order of rxn = number of reactants
    k_stoc = np.copy(k_det)

    for ind in range(nr):
        # If highest order is 3
        if np.max(V_r[ind, :]) == 3:
            k_stoc[ind] = k_det[ind] * 6 / np.power(Na * volume, 2)
        elif np.max(V_r[ind, :]) == 2:  # Highest order is 2
            k_stoc[ind] = k_det[ind] * 2 / np.power(Na * volume, orders[ind]-1)
        else:
            k_stoc[ind] = k_det[ind] / np.power(Na * volume, orders[ind]-1)

    return k_stoc
