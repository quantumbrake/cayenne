"""
    The main class for running stochastic simulation
"""

import multiprocessing as mp
from functools import partial
from typing import List, Optional
from warnings import warn

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

from .algorithms.direct import direct
from .algorithms.tau_adaptive import tau_adaptive
from .algorithms.tau_leaping import tau_leaping
from .model_io import ModelIO
from .results import Results


def wrapper(x, func):
    return func(*x)


class Simulation:
    """
        A main class for running simulations.

        Parameters
        ----------
        species_names : List[str]
            List of species names
        rxn_names : List[str]
            List of reaction names
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
        volume: float, optional
            The volume of the reactor vessel which is important for second
            and higher order reactions. Defaults to 1 arbitrary units.
        chem_flag: bool, optional
            If True, divide by Na (Avogadro's constant) while calculating
            stochastic rate constants. Defaults to ``False``.

        Attributes
        ----------
        results: Results
            The results instance

        Raises
        ------
        ValueError
            If supplied with order > 3.

        Examples
        --------
        >>> V_r = np.array([[1,0],[0,1],[0,0]])
        >>> V_p = np.array([[0,0],[1,0],[0,1]])
        >>> X0 = np.array([10,0,0])
        >>> k = np.array([1,1])
        >>> sim = Simulation(V_r, V_p, X0, k)

        Notes
        -----
        Stochastic reaction rates depend on the size of the system for second
        and third order reactions. By this, we mean the volume of the system in
        which the reactants are contained. Intuitively, this makes sense
        considering that collisions between two or more molecules becomes
        rarer as the size of the system increases. A detailed mathematical
        treatment of this idea can be found in [3]_ .

        In practice, this means that ``volume`` and ``chem_flag`` need to be
        supplied for second and third order reactions. ``volume``
        represents the size of the system containing the reactants.

        In chemical systems ``chem_flag`` should generally be set to ``True``
        as ``k_det`` is specified in units of molarity or M or mol/L.
        For example, a second order rate constant could be = 0.15 mol / (L s).
        Then Avogadro's constant (:math:`N_a`) is used for normalization while
        computing ``k_stoc`` (:math:`c_\mu` in [3]_ ) from ``k_det``.

        In biological systems, ``chem_flag`` should be generally be set to
        ``False`` as ``k_det`` is specified in units of copies/L or CFU/L.
        For example, a second order rate constant could be = 0.15 CFU / (L s).

        References
        ----------
        .. [3] Gillespie, D.T., 1976.
            A general method for numerically simulating the stochastic time
            evolution of coupled chemical reactions. J. Comput. Phys. 22, 403–434.
            doi:10.1016/0021-9991(76)90041-3.

    """

    _results: Optional[Results] = None

    def __init__(
        self,
        species_names: List[str],
        rxn_names: List[str],
        react_stoic: np.ndarray,
        prod_stoic: np.ndarray,
        init_state: np.ndarray,
        k_det: np.ndarray,
        chem_flag: bool = False,
        volume: float = 1.0,
    ) -> None:
        self.species_names = species_names
        self.rxn_names = rxn_names
        self._react_stoic = react_stoic
        self._prod_stoic = prod_stoic
        self._init_state = init_state
        self._k_det = k_det
        self._chem_flag = chem_flag
        self._ns = self._react_stoic.shape[0]
        self._nr = self._react_stoic.shape[1]
        self._volume = volume
        self._orders = np.sum(
            self._react_stoic, axis=0
        )  # Order of rxn = number of reactants
        self._check_consistency()

    def _check_consistency(self):
        if len(self.species_names) != self._ns:
            raise ValueError(
                "Species names must match the stoichiometric matrix shapes."
            )
        if len(self.rxn_names) != self._nr:
            raise ValueError(
                "Reaction names must match the stoichiometric matrix shapes."
            )
        if (self._ns != self._prod_stoic.shape[0]) or (
            self._nr != self._prod_stoic.shape[1]
        ):
            raise ValueError("react_stoic and prod_stoic should be the same shape.")
        if np.any(self._react_stoic < 0):
            raise ValueError("react_stoic cannot have negative elements.")
        if np.any(self._prod_stoic < 0):
            raise ValueError("V_p cannot have negative elements.")
        if np.any(self._init_state < 0):
            raise ValueError("Initial numbers in X0 can't be negative.")
        if np.any(self._k_det < 0):
            raise ValueError("Rate constant(s) can't be negative.")
        if self._k_det.shape[0] != self._nr:
            raise ValueError("Number of rate constants must equal number of reactions")
        if len(self._k_det.shape) != 1:
            raise ValueError("k_det must be a 1-D array")
        if self._chem_flag not in (True, False):
            raise ValueError("chem_flag must be a boolean True or False.")
        if np.max(self._orders) > 3:
            raise ValueError("Order greater than 3 not suppported.")
        if self._ns != self._init_state.shape[0]:
            raise ValueError(
                "X0 must have be of length = num. of species (or rows of V_r)"
            )
        if len(self._init_state.shape) != 1:
            raise ValueError("X0 must be a 1-D array")

    @property
    def results(self) -> Optional[Results]:
        """
            The ``Results`` instance of the simulation

            Returns
            -------
            Optional[Results]
        """
        if self._results is None:
            warn(
                "Run `Simulation.simulate` before requesting the results object",
                Warning,
            )
            return self._results
        else:
            return self._results

    @classmethod
    def load_model(cls, contents: str, contents_type: str) -> "Simulation":
        """
            Load model contents into a Simulation object

            Parameters
            ----------
            model_contents : str
                Either the model string or the file path
            content_type : str, {"ModelString", "ModelFile"}
                The type of the model

            Returns
            -------
            sim : Simulation
                An instance of the Simulation class.
         """
        modelio = ModelIO(contents, contents_type)
        (
            species_names,
            rxn_names,
            react_stoic,
            prod_stoic,
            init_state,
            k_det,
            chem_flag,
            volume,
        ) = modelio.args
        return cls(
            species_names,
            rxn_names,
            react_stoic,
            prod_stoic,
            init_state,
            k_det,
            chem_flag,
            volume,
        )

    def simulate(
        self,
        max_t: float = 10.0,
        max_iter: int = 1000,
        seed: int = 0,
        n_rep: int = 1,
        n_procs: Optional[int] = 1,
        algorithm: str = "direct",
        debug: bool = False,
        **kwargs,
    ) -> None:
        """
            Run the simulation

            Parameters
            ----------
            max_t: float, optional
                The end time of the simulation.
                The default is 10.0.
            max_iter: int, optional
                The maximum number of iterations of the simulation loop.
                The default is 1000 iterations.
            seed: int, optional
                The seed used to generate simulation seeds.
                The default value is 0.
            n_rep: int, optional
                The number of repetitions of the simulation required.
                The default value is 1.
            n_procs: int, optional
                The number of cpu cores to use for the simulation.
                Use ``None`` to automatically detect number of cpu cores.
                The default value is 1.
            algorithm: str, optional
                The algorithm to be used to run the simulation.
                The default value is ``"direct"``.

            Notes
            -----
            The status indicates the status of the simulation at exit. Each
            repetition will have a status associated with it, and these are
            accessible through the ``Simulation.results.status_list``.

            status: int
                Indicates the status of the simulation at exit.

                1: Succesful completion, terminated when ``max_iter`` iterations reached.

                2: Succesful completion, terminated when ``max_t`` crossed.

                3: Succesful completion, terminated when all species went extinct.

                -1: Failure, order greater than 3 detected.

                -2: Failure, propensity zero without extinction.
        """
        tlist = []
        xlist = []
        status_list = []

        if not isinstance(seed, int):
            raise TypeError("Seed should be of type int")
        np.random.seed(seed)
        sim_seeds = np.random.randint(low=0, high=1e7, size=n_rep)

        if not isinstance(max_iter, int):
            raise TypeError("max_iter should be of type int")

        algo_args = []
        if algorithm == "direct":
            for index in range(n_rep):
                algo_args.append(
                    (
                        self._react_stoic,
                        self._prod_stoic,
                        self._init_state,
                        self._k_det,
                        max_t,
                        max_iter,
                        self._volume,
                        sim_seeds[index],
                        self._chem_flag,
                    )
                )
                algo = direct
        elif algorithm == "tau_leaping":
            if "tau" in kwargs.keys():
                tau = kwargs["tau"]
            else:
                tau = 0.1
            for index in range(n_rep):
                algo_args.append(
                    (
                        self._react_stoic,
                        self._prod_stoic,
                        self._init_state,
                        self._k_det,
                        tau,
                        max_t,
                        self._volume,
                        sim_seeds[index],
                        self._chem_flag,
                    )
                )
                algo = tau_leaping
        elif algorithm == "tau_adaptive":
            if "epsilon" in kwargs.keys():
                epsilon = kwargs["epsilon"]
            else:
                epsilon = 0.03
            if "nc" in kwargs.keys():
                nc = kwargs["nc"]
            else:
                nc = 10
            hor = self.HOR
            for index in range(n_rep):
                algo_args.append(
                    (
                        self._react_stoic,
                        self._prod_stoic,
                        self._init_state,
                        self._k_det,
                        hor.astype(np.int),
                        nc,
                        epsilon,
                        max_t,
                        max_iter,
                        self._volume,
                        sim_seeds[index],
                        self._chem_flag,
                    )
                )
            algo = tau_adaptive
        else:
            raise ValueError("Requested algorithm not supported")
        algo_func = partial(wrapper, func=algo)
        if debug:
            results = map(algo_func, algo_args)
        else:
            with mp.Pool(processes=n_procs) as pool:
                results = pool.map(algo_func, algo_args)
        for t, X, status in results:
            tlist.append(t)
            xlist.append(X)
            status_list.append(status)
        self._results = Results(
            self.species_names,
            self.rxn_names,
            tlist,
            xlist,
            status_list,
            algorithm,
            sim_seeds,
        )

    @property
    def HOR(self) -> np.ndarray:
        """
            Determine the HOR vector. HOR(i) is the highest order of reaction
            in which species S_i appears as a reactant.

            Returns
            -------
            HOR: np.ndarray
                Highest order of the reaction for the reactive species as
                defined under Eqn. (27) of [1]_. HOR can be 1, 2 or 3
                if the species appears only once in the reactants.
                If HOR is -2, it appears twice in a second order reaction.
                If HOR is -3, it appears thrice in a third order reaction.
                If HOR is -32, it appears twice in a third order reaction.
                The corresponding value of ``g_i`` in Eqn. (27) is handled
                by ``tau_adaptive``.

            References
            ----------
            .. [1] Cao, Y., Gillespie, D.T., Petzold, L.R., 2006.
                Efficient step size selection for the tau-leaping simulation
                method. J. Chem. Phys. 124, 044109. doi:10.1063/1.2159468
        """
        ns = self._react_stoic.shape[0]
        HOR = np.zeros(ns, dtype=np.int32)
        orders = np.sum(self._react_stoic, axis=0)
        for ind in range(ns):
            this_orders = orders[np.where(self._react_stoic[ind, :] > 0)[0]]
            if len(this_orders) == 0:
                HOR[ind] = 0
                continue
            HOR[ind] = np.max(this_orders)
            if HOR[ind] == 1:
                continue
            order_2_indices = np.where(orders == 2)
            this_react_stoic = self._react_stoic[ind, :]
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

    def plot(self, species_names: list = None, new_names: list = None):
        """
            Plot the simulation

            Parameters
            ----------
            species_names: list, optional
                The names of the species to be plotted (``list`` of ``str``).
                The default is ``None`` and  plots all species.
            new_names: list, optional
                The names of the species to be plotted.
                The default is ``"xi"`` for species ``i``.

            Returns
            -------
            fig: class 'matplotlib.figure.Figure'
                Figure object of the generated plot.
            ax: class 'matplotlib.axes._subplots.AxesSubplot
                Axis objected of the generated plot.
        """

        if self._results is None:
            raise ValueError("Simulate not run.")
        else:
            if species_names is None:
                species_names = self.species_names
            n_species = len(species_names)
            prop_cycle = plt.rcParams["axes.prop_cycle"]
            colors = prop_cycle.by_key()["color"]
            fig, ax = plt.subplots()
            res = self._results
            legend_handlers = [0] * n_species
            generic_names = [""] * n_species
            for index1 in range(n_species):
                legend_handlers[index1] = mlines.Line2D([], [], color=colors[index1])
                this_species = species_names[index1]
                generic_names[index1] = this_species
                n_reps = len(res)
                for index2 in range(n_reps):
                    ax.step(
                        res.t_list[index2],
                        res.get_species([this_species])[index2],
                        color=colors[index1],
                        where="post",
                    )
            if new_names is None:
                new_names = self.species_names
            fig.legend(legend_handlers, new_names)
            ax.set_xlabel("Time")
            ax.set_ylabel("Species amounts")
            return fig, ax
