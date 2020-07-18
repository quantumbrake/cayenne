"""
    Implementation of Gillespie's Direct method. This is an exact simulation
    algorithm that simulates each reaction step. This makes it slower than
    other methods, but it's a good place to start.
"""

cimport cython
cimport numpy as np
import numpy as np
# import random # faster than using np.random, but reproducibility issues
from ..utils cimport roulette_selection
from ..utils import get_kstoc
from libc.math cimport log


@cython.boundscheck(False)
@cython.wraparound(False)
def direct(
    react_stoic: np.ndarray,
    prod_stoic: np.ndarray,
    init_state: np.ndarray,
    k_det: np.ndarray,
    max_t: float,
    max_iter: int,
    volume: float,
    seed: int,
    chem_flag: bool,
):
    """
        Runs the Direct Stochastic Simulation Algorithm [#]_

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
        max_t: float
            The maximum simulation time to run the simulation for.
        max_iter: int
            The maximum number of iterations to run the simulation for.
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
            Numpy array of the states of the system at times in in ``t``.
        status: int
            Indicates the status of the simulation at exit.

            1 - Succesful completion, terminated when ``max_iter`` iterations reached.

            2 - Succesful completion, terminated when ``max_t`` crossed.

            3 - Succesful completion, terminated when all species went extinct.

            -1 - Failure, order greater than 3 detected.

            -2 - Failure, propensity zero without extinction.

        References
        ----------
        .. [#] Gillespie, D.T., 1976.
            A general method for numerically simulating the stochastic time
            evolution of coupled chemical reactions. J. Comput. Phys. 22, 403–434.
            doi:10.1016/0021-9991(76)90041-3.
    """

    cdef:
        int ite=1, ind=0, ind1=0, ind2=0, continue_flag=0, choice, status
        double t_curr=0, prop_sum=0
        Py_ssize_t ns=react_stoic.shape[0], nr=react_stoic.shape[1]
    v = prod_stoic - react_stoic  # ns x nr
    x = np.zeros((max_iter, ns), dtype=np.int64)
    t = np.zeros((max_iter))
    x[0, :] = init_state.copy()
    xtemp = init_state.copy()  # Temporary X for updating
    status = 0
    np.random.seed(seed)
    # random.seed(seed)  # Set the seed
    # Determine kstoc from kdet and the highest order or reactions
    prop = get_kstoc(react_stoic, k_det, volume, chem_flag)  # Vector of propensities
    kstoc = prop.copy()  # Stochastic rate constants
    cdef double [:] kstoc_view = kstoc
    cdef double [:] prop_view = prop
    cdef long long [:] xtemp_view = xtemp
    cdef long [:, :] v_view = v
    cdef long [:, :] react_stoic_view = react_stoic
    cdef long long [:, :] x_view = x
    while ite < max_iter:
        # Calculate propensities
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
        # Roulette wheel
        roulette_results = roulette_selection(prop_view, x_view[ite-1, :])
        choice = roulette_results[0]
        status = roulette_results[1]
        if status == 0:
            for ind1 in range(ns):
                xtemp_view[ind1] = x_view[ite-1, ind1] + v_view[ind1, choice]
        else:
            return t[:ite], x[:ite, :], status

        # If negative species produced, reject step
        # if np.min(xtemp) < 0:
        #     continue
        for ind in range(ns):
            if xtemp_view[ind] < 0:
                continue_flag = 1
                break
        if continue_flag:
            continue_flag = 0
            continue
        # Update xt and t_curr
        else:
            r2 = np.random.random()
            prop_sum = 0
            for ind in range(nr):
                prop_sum += prop_view[ind]
            # t_curr += 1 / prop_sum * log(1 / r2)
            # t_curr += -1 / prop_sum * log(r2)
            t_curr += np.random.exponential(1/prop_sum)
            if t_curr > max_t:
                status = 2
                return t[:ite], x[:ite, :], status
        for ind in range(nr):
            prop_view[ind] = kstoc_view[ind]
        for ind in range(ns):
            x_view[ite, ind] = xtemp_view[ind]
        t[ite] = t_curr
        ite += 1
    status = 1
    return t[:ite], x[:ite, :], status
