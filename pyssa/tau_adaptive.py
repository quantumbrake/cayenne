"""
    Implementation of the tau leaping algorithm in Numba
"""

from typing import Tuple
from numba import njit, jit
import numpy as np
from .utils import get_kstoc, roulette_selection, HIGH, TINY
from .direct import direct


@njit(nogil=True, cache=False)
def get_HOR(react_stoic: np.ndarray):
    """ Determine the HOR vector. HOR(i) is the highest order of reaction
        in which species S_i appears as a reactant.

        Parameters
        ----------
        react_stoic : (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are rows and species are columns.

        Returns
        -------
        HOR : np.ndarray
            Highest order of the reaction for the reactive species as
            defined under Eqn. (27) of [1]_. HOR can be 1, 2 or 3
            if the species appears only once in the reactants.
            If HOR is -2, it appears twice in a second order reaction.
            If HOR is -3, it appears thrice in a third order reaction.
            If HOR is -32, it appears twice in a third order reaction.
            The corresponding value of `g_i` in Eqn. (27) is handled
            by `tau_adaptive`.

        References
        ----------
        .. [1] Cao, Y., Gillespie, D.T., Petzold, L.R., 2006.
        Efficient step size selection for the tau-leaping simulation
        method. J. Chem. Phys. 124, 044109. doi:10.1063/1.2159468
    """
    ns = react_stoic.shape[0]
    HOR = np.zeros(ns, dtype=np.int64)
    orders = np.sum(react_stoic, axis=0)
    for ind in range(ns):
        this_orders = orders[np.where(react_stoic[ind, :] > 0)[0]]
        if len(this_orders) == 0:
            HOR[ind] = 0
            continue
        HOR[ind] = np.max(this_orders)
        if HOR[ind] == 1:
            continue
        order_2_indices = np.where(orders == 2)
        this_react_stoic = react_stoic[ind, :]
        if order_2_indices[0].size > 0:
            if np.max(this_react_stoic[order_2_indices[0]]) == 2 and HOR[ind] == 2:
                HOR[ind] = -2  # g_i should be (2 + 1/(x_i-1))
        if np.where(orders == 3):
            if (
                HOR[ind] == 3
                and np.max(this_react_stoic[np.where(this_orders == 3)[0]]) == 2
            ):
                HOR[ind] = -32  # g_i should be (3/2 * (2 + 1/(x_i-1)))
            elif (
                HOR[ind] == 3
                and np.max(this_react_stoic[np.where(this_orders == 3)[0]]) == 3
            ):
                HOR[ind] = -3  # g_i should be(3 + 1/(x_i-1) + 2/(x_i-2))
    return HOR


@njit(nogil=True, cache=False)
def step1(kstoc, xt, react_stoic, v, nc):
    """ Determine critical reactions """
    nr = react_stoic.shape[1]
    ns = react_stoic.shape[0]
    prop = np.copy(kstoc)
    L = np.zeros(nr, dtype=np.int64)
    # Calculate the propensities
    for ind1 in range(nr):
        for ind2 in range(ns):
            # prop = kstoc * product of (number raised to order)
            prop[ind1] *= np.power(xt[ind2], react_stoic[ind2, ind1])
    for ind in range(nr):
        vis = v[:, ind]
        L[ind] = np.nanmin(xt[vis < 0] / np.abs(vis[vis < 0]))
    # A reaction j is critical if Lj <nc. However criticality is
    # considered only for reactions with propensity greater than
    # 0 (`prop > 0`).
    crit = (L < nc) * (prop > 0)
    # To get the non-critical reactions, we use the bitwise not operator.
    not_crit = ~crit
    return prop, crit, not_crit


@njit(nogil=True, cache=False)
def step2(not_crit, react_species, v, xt, HOR, prop, epsilon):
    """ 2. Generate candidate taup """
    n_react_species = react_species.shape[0]
    mup = np.zeros(n_react_species, dtype=np.float64)
    sigp = np.zeros(n_react_species, dtype=np.float64)
    tau_num = np.zeros(n_react_species, dtype=np.float64)
    if np.sum(not_crit) == 0:
        taup = HIGH
    else:
        # Compute mu from eqn 32a and sig from eqn 32b
        for ind, species_index in enumerate(react_species):
            this_v = v[species_index, :]
            for i in range(len(not_crit)):
                if not_crit[i]:
                    mup[ind] += this_v[i] * prop[i]
                    sigp[ind] += this_v[i] * prop[i] * this_v[i]
            if mup[ind] == 0:
                mup[ind] = TINY
            if sigp[ind] == 0:
                sigp[ind] = TINY
            if HOR[species_index] > 0:
                g = HOR[species_index]
            elif HOR[species_index] == -2:
                if xt[species_index] is not 1:
                    g = 1 + 2 / (xt[species_index] - 1)
                else:
                    g = 2
            elif HOR[species_index] == -3:
                if xt[species_index] not in [1, 2]:
                    g = 3 + 1 / (xt[species_index] - 1) + 2 / (xt[species_index] - 2)
                else:
                    g = 3
            elif HOR[species_index] == -32:
                if xt[species_index] is not 1:
                    g = 3 / 2 * (2 + 1 / (xt[species_index] - 1))
                else:
                    g = 3
            tau_num[ind] = max(epsilon * xt[species_index] / g, 1)
        taup = np.nanmin(
            np.concatenate((tau_num / np.abs(mup), np.power(tau_num, 2) / np.abs(sigp)))
        )
    return taup


@njit(nogil=True, cache=False)
def step5(taup, taupp, nr, not_crit, prop, xt):
    K = np.zeros(nr, dtype=np.int64)
    if taup < taupp:
        tau = taup
        for ind in range(nr):
            if not_crit[ind]:
                K[ind] = np.random.poisson(prop[ind] * tau)
            else:
                K[ind] = 0
    else:
        tau = taupp
        # Identify the only critical reaction to fire
        # Send in xt to match signature of roulette_selection
        temp = prop.copy()
        temp[not_crit] = 0
        j_crit, _ = roulette_selection(temp, xt)
        for ind in range(nr):
            if not_crit[ind]:
                K[ind] = np.random.poisson(prop[ind] * tau)
            elif ind == j_crit:
                K[ind] = 1
            else:
                K[ind] = 0
    return tau, K


@njit(nogil=True, cache=False)
def tau_adaptive(
    react_stoic: np.ndarray,
    prod_stoic: np.ndarray,
    init_state: np.ndarray,
    k_det: np.ndarray,
    nc: int,
    epsilon: float,
    max_t: float,
    max_iter: int,
    volume: float,
    seed: int,
    chem_flag: bool,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
        Parameters
        ---------
        react_stoic : (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are columns and species are rows.
        prod_stoic : (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the products.
            Reactions are columns and species are rows.
        init_state : (ns,) ndarray
            A 1D array representing the initial state of the system.
        k_det : (nr,) ndarray
            A 1D array representing the deterministic rate constants of the
            system.
        nc : int
            The criticality threshold. Reactions with that cannot fire more than
            nc times are deemed critical.
        epsilon : float
            The epsilon used in tau-leaping, measure of the bound on relative
            change in propensity.
        max_t : float
            The maximum simulation time to run the simulation for.
        max_iter : int
            The maximum number of iterations to run the simulation for.
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
    ns = react_stoic.shape[0]
    nr = react_stoic.shape[1]
    v = prod_stoic - react_stoic  # ns x nr
    x = np.zeros((max_iter, ns), dtype=np.int64)
    xt = np.zeros(ns, dtype=np.int64)
    t = np.zeros((max_iter))
    x[0, :] = init_state.copy()
    np.random.seed(seed)  # Set the seed
    # Determine kstoc from kdet and the highest order or reactions
    prop = np.copy(
        get_kstoc(react_stoic, k_det, volume, chem_flag)
    )  # Vector of propensities
    kstoc = prop.copy()  # Stochastic rate constants
    prop_sum = np.sum(prop)

    # Determine the HOR vector. HOR(i) is the highest order of reaction
    # in which species S_i appears as a reactant.
    HOR = get_HOR(react_stoic)

    if np.sum(prop) < TINY:
        if np.sum(x[ite - 1, :]) > TINY:
            status = -2
            return t[:ite], x[:ite, :], status

    M = nr
    N = ns
    vis = np.zeros(M, dtype=np.int64)
    react_species = np.where(np.sum(react_stoic, axis=1) > 0)[0]
    skip_flag = False

    while ite < max_iter:

        # If negatives are not detected in step 6, perform steps 1 and 2.
        # Else skip to step 3.
        if not skip_flag:
            xt = x[ite - 1, :]

            # Step 1:
            prop, crit, not_crit = step1(kstoc, xt, react_stoic, v, nc)
            prop_sum = np.sum(prop)
            if prop_sum < TINY:
                status = 3
                return t[:ite], x[:ite, :], status

            # Step 2:
            taup = step2(not_crit, react_species, v, xt, HOR, prop, epsilon)

        # 3. For small taup, do SSA
        # -------------------------
        skip_flag = False
        if taup < 10 / prop_sum:
            t_ssa, x_ssa, status = direct(
                react_stoic,
                prod_stoic,
                x[ite - 1, :],
                k_det,
                max_t=max_t - t[ite - 1],
                max_iter=min(101, max_iter - ite),
                volume=volume,
                seed=seed,
                chem_flag=chem_flag,
            )
            # t_ssa first element is 0. x_ssa first element is x[ite - 1, :].
            # Both should be dropped while logging the results.
            len_simulation = len(t_ssa) - 1  # Since t_ssa[0] is 0
            t[ite : ite + len_simulation] = t_ssa[1:] + t[ite-1]
            x[ite : ite + len_simulation, :] = x_ssa[1:]
            ite += len_simulation
            if status == 3 or status == 2:
                return t, x, status
            continue

        # 4. Generate second candidate taupp
        # ----------------------------------
        if np.sum(prop[crit]) == 0:
            taupp = HIGH
        else:
            taupp = 1 / np.sum(prop[crit]) * np.log(1 / np.random.rand())

        # 5. Leap
        # -------
        tau, K = step5(taup, taupp, nr, not_crit, prop, xt)

        # 6. Handle negatives, update and exit conditions
        # -----------------------------------------------
        # v = ns, nr
        # K = nr
        vdotK = np.zeros((ns,), dtype=np.int64)
        for ind1 in range(ns):
            for ind2 in range(nr):
                vdotK[ind1] += v[ind1, ind2] * K[ind2]
        x_new = x[ite - 1, :] + vdotK
        if np.any(x_new < 0):
            taup = taup / 2
            skip_flag = True
        else:
            # Update states if nothing is negative
            x[ite, :] = x[ite - 1, :] + vdotK
            t[ite] = t[ite - 1] + tau
            ite += 1

        # Exit conditions
        if t[ite - 1] > max_t:
            status = 2
            return t[:ite], x[:ite, :], status

    status = 1
    return t[:ite], x[:ite, :], status
