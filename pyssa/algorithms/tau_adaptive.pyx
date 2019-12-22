"""
    Experimental implementation of the `tau adaptive algorithm
    <https://doi.org/10.1063/1.2159468>`_. This is an approximate method that
    builds off the tau_leaping method. It self-adapts the value of ``tau``
    over the course of the simulation. For sytems with a small number of
    molecules, it will be similar in speed to the direct method. For systems
    with a large number of molecules, it will be much faster than the direct
    method.
"""

cimport cython
cimport numpy as np
import numpy as np
# import random
from ..utils import get_kstoc, TINY, HIGH
from ..utils cimport roulette_selection
from .direct import direct
from libc.math cimport log


cdef step1(
        double [:] kstoc_view,
        long long [:] xt_view,
        long [:, :] react_stoic_view,
        long [:, :] v_view,
        int nc,
    ):
    """ Determine critical reactions """
    cdef:
        Py_ssize_t ns=react_stoic_view.shape[0], nr=react_stoic_view.shape[1]
        double prop_sum = 0
        long [:] vis
        int visflag = 0
    cdef double [:] prop_view = np.zeros(nr)
    prop_view[...] = kstoc_view
    L = np.ones(nr, dtype=np.int64) * (nc + 1)
    cdef np.int32_t [:] crit = np.zeros(nr, dtype=np.int32)
    cdef np.int32_t [:] not_crit = np.zeros(nr, dtype=np.int32)
    # Calculate the propensities
    for ind1 in range(nr):
        for ind2 in range(ns):
            # prop = kstoc * product of (number raised to order)
            if react_stoic_view[ind2, ind1] == 1:
                prop_view[ind1] *= xt_view[ind2]
            elif react_stoic_view[ind2, ind1] == 2:
                prop_view[ind1] *= xt_view[ind2] * (xt_view[ind2] - 1) / 2
            elif react_stoic_view[ind2, ind1] == 3:
                prop_view[ind1] *= xt_view[ind2] * (xt_view[ind2] - 1) * (xt_view[ind2] - 2) / 6
    for ind1 in range(nr):
        vis = v_view[:, ind1]
        for ind2 in range(ns):
            if vis[ind2] < 0:
                L[ind1] = min(L[ind1], xt_view[ind2] / abs(vis[ind2]))
    # A reaction j is critical if Lj <nc. However criticality is
    # considered only for reactions with propensity greater than
    # 0 (`prop > 0`).
    for ind1 in range(nr):
        if (L[ind1] < nc) and (prop_view[ind1] > 0):
            crit[ind1] = 1
            not_crit[ind1] = 0
        else:
            crit[ind1] = 0
            not_crit[ind1] = 1
    # crit = (L < nc) * (prop > 0)
    return prop_view, crit, not_crit


def _py_step1(kstoc, xt, react_stoic, v, nc):
    prop_view, crit, not_crit = step1(kstoc, xt, react_stoic, v, nc)
    return np.array(prop_view), np.array(crit), np.array(not_crit)


cdef inline step2_get_g(int hor, long long x):
    cdef double g
    if hor > 0:
        g = float(hor)
    elif hor == -2:
        if x is not 1:
            g = 2.0 + 1.0 / (x - 1.0)
        else:
            g = 2.0
    elif hor == -3:
        if x not in [1, 2]:
            g = 3.0 + 1.0 / (x - 1.0) + 2.0 / (x - 2.0)
        else:
            g = 3.0
    elif hor == -32:
        if x is not 1:
            g = 3.0 / 2.0 * (2.0 + 1.0 / (x - 1.0))
        else:
            g = 3.0
    return g


def _py_step2_get_g(hor, x):
    return step2_get_g(hor, x)

cdef step2(
    int [:] not_crit,
    np.int32_t [:] react_species_view,
    long [:, :] v_view,
    long long [:] xt_view,
    np.int32_t [:] hor_view,
    double [:] prop,
    double epsilon,
    ):

    """ 2. Generate candidate taup """
    cdef:
        Py_ssize_t n_react_species = react_species_view.shape[0]
        Py_ssize_t ns=v_view.shape[0], nr=v_view.shape[1]
    cdef double [:] mup_view = np.zeros(n_react_species, dtype=np.float64)
    cdef double [:] sigp_view = np.zeros(n_react_species, dtype=np.float64)
    cdef double [:] tau_num = np.zeros(n_react_species, dtype=np.float64)
    cdef double taup = HIGH
    cdef int crit_flag = 0
    for ind in range(nr):
        if not_crit[ind] != 0:
            crit_flag = 1
            break
    if not crit_flag:
        taup = HIGH
    else:
        # Compute mu from eqn 32a and sig from eqn 32b
        for ind, species_index in enumerate(react_species_view):
            for i in range(nr):
                if not_crit[i]:
                    mup_view[ind] += v_view[species_index, i] * prop[i]
                    sigp_view[ind] += v_view[species_index, i] * prop[i] * v_view[species_index, i]
            if mup_view[ind] == 0:
                mup_view[ind] = TINY
            if sigp_view[ind] == 0:
                sigp_view[ind] = TINY
            g = step2_get_g(hor_view[species_index], xt_view[species_index])
            tau_num[ind] = max(epsilon * xt_view[species_index] / g, 1)
        for ind in range(n_react_species):
            if mup_view[ind] != 0:
                v1 = tau_num[ind] / abs(mup_view[ind])
            else:
                v1 = HIGH
            if sigp_view[ind] != 0:
                v2 = tau_num[ind] * tau_num[ind] / sigp_view[ind]
            else:
                v2 = HIGH
            taup = min(taup, min(v1, v2))
        # print("mup, sigp, tau_num, xt, taup", np.array(mup_view), np.array(sigp_view), np.array(tau_num), np.array(xt_view), taup)
    return taup


def _py_step2(not_crit,
    react_species,
    v,
    xt,
    hor,
    prop,
    epsilon):
    taup = step2(not_crit, react_species, v, xt, hor, prop, epsilon)
    return np.array(taup)


cdef step5(
    double taup,
    double taupp,
    Py_ssize_t nr,
    int [:] not_crit,
    double [:] prop_view,
    long long [:] xt_view,
    ):
    cdef:
        double tau
        double [:] temp_prop = np.zeros(nr)
        long long [:] K = np.zeros(nr, dtype=np.int64)
    if taup < taupp:
        tau = taup
        for ind in range(nr):
            if not_crit[ind]:
                K[ind] = np.random.poisson(prop_view[ind] * tau)
            else:
                K[ind] = 0
        # print("taup, prop_view, K", taup, np.array(prop_view), np.array(K))
    else:
        tau = taupp
        # Identify the only critical reaction to fire
        # Send in xt to match signature of roulette_selection
        for ind in range(nr):
            if not_crit[ind]:
                temp_prop[ind] = 0.0
            else:
                temp_prop[ind] = prop_view[ind]
        j_crit, _ = roulette_selection(temp_prop, xt_view)
        for ind in range(nr):
            if not_crit[ind]:
                K[ind] = np.random.poisson(prop_view[ind] * tau)
            elif ind == j_crit:
                K[ind] = 1
            else:
                K[ind] = 0
    return tau, K


def _py_step5(taup, taupp, nr, not_crit, prop, xt):
    tau, K = step5(taup, taupp, nr, not_crit, prop, xt)
    return tau, np.array(K)

@cython.boundscheck(False)
@cython.wraparound(False)
def tau_adaptive(
    react_stoic: np.ndarray,
    prod_stoic: np.ndarray,
    init_state: np.ndarray,
    k_det: np.ndarray,
    hor: np.ndarray,
    nc: int,
    epsilon: float,
    max_t: float,
    max_iter: int,
    volume: float,
    seed: int,
    chem_flag: bool,
):
    """
        Runs the adaptive tau leaping simulation algorithm [#]_

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
        hor
            A 1D array of the highest order reaction in which each species
            appears.
        nc: int
            The criticality threshold. Reactions with that cannot fire more than
            nc times are deemed critical.
        epsilon: float
            The epsilon used in tau-leaping, measure of the bound on relative
            change in propensity.
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

            -3 - Negative species count encountered

        References
        ----------
        .. [#] Cao, Y., Gillespie, D.T., Petzold, L.R., 2006.
            Efficient step size selection for the tau-leaping simulation
            method. J. Chem. Phys. 124, 044109. doi:10.1063/1.2159468
    """
    cdef:
        int ite = 1
        long long xtsum = 0
        bint skipflag = 0, negflag = 0
        Py_ssize_t ns=react_stoic.shape[0], nr=react_stoic.shape[1]
        double prop_sum = 0, prop_crit_sum = 0, taup, taupp
    v = prod_stoic - react_stoic  # ns x nr
    x = np.zeros((max_iter, ns), dtype=np.int64)
    t = np.zeros((max_iter))
    x[0, :] = init_state.copy()
    # random.seed(seed)  # Set the seed
    np.random.seed(seed)
    # Determine kstoc from kdet and the highest order or reactions
    prop = get_kstoc(react_stoic, k_det, volume, chem_flag) # Vector of propensities
    kstoc = prop.copy()  # Stochastic rate constants
    vis = np.zeros(nr, dtype=np.int64)
    react_species = np.where(np.sum(react_stoic, axis=1) > 0)[0]
    react_species = react_species.astype(dtype=np.int32)

    cdef double [:] kstoc_view = kstoc
    cdef double [:] prop_view = prop
    cdef long [:, :] v_view = v
    cdef long [:, :] react_stoic_view = react_stoic
    cdef long long [:, :] x_view = x
    cdef double [:] t_view = t
    cdef np.int32_t [:] react_species_view = react_species
    cdef np.int32_t [:] hor_view = hor.astype(np.int32)
    cdef long long [:] vdotK = np.zeros((ns,), dtype=np.int64)
    cdef long long [:] x_new = np.zeros((ns,), dtype=np.int64)


    for ind in range(nr):
        prop_sum += prop_view[ind]
    if prop_sum < TINY:
        for ind in range(ns):
            xtsum += x_view[0, ind]
        if xtsum > TINY:
            status = -2
            return t[:ite], x[:ite, :], status


    while ite < max_iter:
        if not skipflag:

            # Step 1
            # print("starting step 1", ite, t[ite-1])
            prop_view, crit, not_crit = step1(
                kstoc_view,
                x_view[ite - 1, :],
                react_stoic_view,
                v_view,
                nc
            )
            prop_sum = 0
            for ind in range(nr):
                prop_sum += prop_view[ind]
            if prop_sum < TINY:
                status = 3
                return t[:ite], x[:ite, :], status

            # Step 2
            # print("starting step 2", ite, t[ite-1])
            taup = step2(
                not_crit,
                react_species_view,
                v_view,
                x_view[ite-1, :],
                hor_view,
                prop_view,
                epsilon,
            )
            # print("ite, taup, x, 10/propsum, time is: ", ite, taup, np.array(x_view[ite-1, :]), 10.0/prop_sum, t[ite-1])

        # Step 3: For small taup, do SSA
        skipflag = 0
        if taup < 10.0 / prop_sum:
            t_ssa, x_ssa, status = direct(
                react_stoic,
                prod_stoic,
                x[ite - 1, :],
                k_det,
                max_t=max_t - t[ite - 1],
                max_iter=min(101, max_iter - ite),
                volume=volume,
                seed=np.random.randint(low=0, high=1e7),
                chem_flag=chem_flag,
            )
            # t_ssa first element is 0. x_ssa first element is x[ite - 1, :].
            # Both should be dropped while logging the results.
            len_simulation = len(t_ssa) - 1  # Since t_ssa[0] is 0
            for ind1 in range(len_simulation):
                t_view[ite+ind1] = t_ssa[1+ind1] + t[ite-1]
                for ind2 in range(ns):
                    x_view[ite+ind1, ind2] = x_ssa[1+ind1, ind2]
            ite += len_simulation
            if status == 3 or status == 2:
                return t[:ite], x[:ite, :], status
            continue

        # Step 4. Generate second candidate taupp
        # ----------------------------------
        # print("starting step 4", ite, t[ite-1])
        prop_crit_sum = 0.0
        for ind in range(nr):
            if crit[ind]:
                prop_crit_sum += prop[crit[ind]]
        if prop_crit_sum == 0:
            taupp = HIGH
        else:
            # taupp = 1 / prop_crit_sum * log(1 / random.rand())
            taupp = 1 / prop_crit_sum * log(1 / np.random.random())
            # taupp = np.random.exponential(1/prop_crit_sum)

        # Step 5. Leap
        # ------------
        tau, K = step5(taup, taupp, nr, not_crit, prop_view, x_view[ite-1, :])
        # print("taup, taupp, tau", taup, taupp, tau)

        # 6. Handle negatives, update and exit conditions
        # -----------------------------------------------
        # v = ns, nr
        # K = nr
        for ind1 in range(ns):
                vdotK[ind1] = 0
        negflag = 0
        for ind1 in range(ns):
            for ind2 in range(nr):
                vdotK[ind1] += v[ind1, ind2] * K[ind2]
            x_new[ind1] = x[ite - 1, ind1] + vdotK[ind1]
            if x_new[ind1] < 0:
                negflag = 1
        # print("vdotk, k, negflag", np.array(vdotK), np.array(K), negflag)
        if negflag:
            taup = taup / 2
            skip_flag = True
        else:
            # Update states if nothing is negative
            for ind1 in range(ns):
                x_view[ite, ind1] = x_view[ite - 1, ind1] + vdotK[ind1]
            t_view[ite] = t_view[ite - 1] + tau
            ite = ite + 1

        # Exit conditions
        if t[ite - 1] > max_t:
            status = 2
            return t[:ite], x[:ite, :], status

    status = 1
    return t[:ite], x[:ite, :], status
