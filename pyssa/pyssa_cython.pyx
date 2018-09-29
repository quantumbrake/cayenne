"""
    Naive implementation of the Gillespie algorithm (direct method) in Cython
"""

import cython
import numpy as np


Na = 6.023e23  # Avogadro's constant


@cython.boundscheck(False)
@cython.wraparound(False)
cdef (int, int) cy_roulette_selection(double[:] prop_list, long[:] Xt):
    """Perform roulette selection on the list of propensities"""
    cdef int prop0 = np.sum(prop_list)  # Sum of propensities
    cdef int status
    if prop0 == 0:
        if np.sum(Xt) == 0:
            status = 3
            return -1, status
        else:
            status = -2
            return -1, status
    cdef double[:] prop = np.divide(prop_list, prop0)  # Normalize propensities to be < 1
    # Concatenate 0 to list of probabilities
    probs = [0] + list(np.cumsum(prop))
    cdef float r1 = np.random.random() # Roll the wheel
    # Identify where it lands and update that reaction
    cdef int ind1
    cdef int choice = 0
    for ind1 in range(len(probs)):
        if r1 <= probs[ind1]:
            choice = ind1 - 1
            break
    return choice, 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:] cy_get_kstoc(double[:] k_det, long[:, :] V_r, float volume = 1.0):
    """Compute k_stoc from k_det"""
    cdef double Na = 6.023e23  # Avogadro's constant
    cdef int nr = V_r.shape[0]  # Number of reactions
    cdef long[:] orders = np.sum(V_r, 1)  # Order of rxn = number of reactants
    k_stoc = np.zeros_like(k_det)
    cdef int ind
    for ind in range(nr):
        # If highest order is 3
        if np.max(V_r[ind, :]) == 3:
            k_stoc[ind] = k_det[ind] * 6 / np.power(Na * volume, 2)
        elif np.max(V_r[ind, :]) == 2:  # Highest order is 2
            k_stoc[ind] = k_det[ind] * 2 / np.power(Na * volume, orders[ind] - 1)
        else:
            k_stoc[ind] = k_det[ind] / np.power(Na * volume, orders[ind] - 1)
    return k_stoc


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cy_direct_naive(
    long[:, :] V_r,
    long[:, :] V_p,
    long[:] X0,
    double[:] k_det,
    float max_t = 1.0,
    long max_iter = 100,
    float volume = 1.0,
    int seed = 0
):
    """Naive implementation of the Direct method"""
    cdef int ite = 1  # Iteration counter
    cdef double t = 0.0  # Time in seconds
    cdef int nr = V_r.shape[0]  # Number of reactions
    cdef int ns = V_r.shape[1]  # Number of species

    if (nr != V_p.shape[0]) or (ns != V_p.shape[1]):
        raise ValueError('V_r and V_p should be the same shape.')

    if (nr != k_det.shape[0]):
        raise ValueError('Number of elements in k_det must equal\
         number of rows in V_r.')

    if np.any(np.less(V_r, 0)):
        raise ValueError('V_r cannot have negative elements.')

    if np.any(np.less(V_p, 0)):
        raise ValueError('V_p cannot have negative elements.')

    if np.any(np.less(X0, 0)):
        raise ValueError('Initial numbers in X0 can\'t be negative.')

    if np.any(np.less(k_det, 0)):
        neg_indices = np.where(np.less(k_det, 0))[0]
        raise ValueError('Rate constant(s) at position(s) ' + str(neg_indices) + ' are negative.')

    V = np.subtract(V_p, V_r)  # nr x ns
    cdef long[:] Xt = np.copy(X0)  # Number of species at time t
    cdef long[:] Xtemp = np.zeros(nr, dtype=int)  # Temporary X for updating
    cdef double[:] kstoc = np.zeros(nr)  # Stochastic rate constants
    cdef long[:] orders = np.sum(V_r, 1)  # Order of rxn = number of reactants
    cdef int status = 0
    np.random.seed(seed=seed)  # Set the seed

    if np.max(orders) > 3:
        raise ValueError('Order greater than 3 detected.')

    # Determine kstoc from kdet and the highest order or reactions
    kstoc = cy_get_kstoc(k_det, V_r, volume)
    prop = np.copy(kstoc)  # Vector of propensities

    while ite < max_iter:
        # Calculate propensities
        for ind1 in range(nr):
            for ind2 in range(ns):
                # prop = kstoc * product of (number raised to order)
                prop[ind1] *= np.power(Xt[ind2], V_r[ind1, ind2])
        # Roulette wheel
        [choice, status] = cy_roulette_selection(prop, Xt)
        if status == 0:
            Xtemp = Xt + V[choice, :]
        else:
            return t, Xt, status

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
                print("Reached maximum time (t = )", t)
                return t, Xt, status
        prop = np.copy(kstoc)
        ite += 1
    status = 1
    return t, Xt, status
