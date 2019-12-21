"""
    Module that defines the `Results` class
"""

from collections.abc import Collection
from typing import List, Tuple
from warnings import warn

import numpy as np


class Results(Collection):
    """
        A class that stores simulation results and provides methods to access them

        Parameters
        ----------
        t_list: List[float]
            List of time points for each repetition
        x_list: List[np.ndarray]
            List of system states for each repetition
        status_list: List[int]
            List of return status for each repetition
        algorithm: str
            Algorithm used to run the simulation
        sim_seeds: List[int]
            List of seeds used for the simulation
    """

    def __init__(
        self,
        t_list: List[np.ndarray],
        x_list: List[np.ndarray],
        status_list: List[int],
        algorithm: str,
        sim_seeds: List[int],
    ) -> None:
        self.x_list = x_list
        self.t_list = t_list
        self.status_list = status_list
        self.algorithm = algorithm
        self.sim_seeds = sim_seeds
        if not self._check_consistency():
            raise ValueError("Inconsistent results passed")

    def _check_consistency(self) -> bool:
        """
            Check consistency of results

            Returns
            -------
            bool
                True if results are consistent
                False otherwise
        """
        if (
            len(self.x_list)
            == len(self.t_list)
            == len(self.status_list)
            == len(self.sim_seeds)
        ):
            pass
        else:
            return False
        for x, t, status in self:
            if x.shape[0] != t.shape[0]:
                return False
            if not isinstance(status, int):
                return False
        return True

    def __repr__(self) -> str:
        """
            Return summary of simulation.

            Returns
            -------
            summary: str
                Summary of the simulation with length of simulation, algorithm and seeds used.
        """
        summary = f"<Results n_rep={len(self)} algorithm={self.algorithm} sim_seeds={self.sim_seeds}>"
        return summary

    def __str__(self) -> str:
        """ Return self.__repr__() """
        return self.__repr__()

    def __iter__(self):
        """ Iterate over each repetition """
        return zip(self.x_list, self.t_list, self.status_list)

    def __len__(self):
        """
            Return number of repetitions in simulation

            Returns
            -------
            n_rep: int
                Number of repetitions in simulation
        """
        n_rep = len(self.x_list)
        return n_rep

    def __contains__(self, ind):
        """ Returns True if ind is one of the repetition numbers """
        if ind < len(self):
            return True
        else:
            return False

    def __getitem__(self, ind: int):
        """
            Return sim. state, time points and status of repetition no. `ind`

            Parameters
            ----------
            ind: int
                Index of the repetition in the simulation

            Returns
            -------
            x_ind: np.ndarray
                Simulation status of repetition no. `ind`
            t_ind: np.ndarray
                Time points of repetition no. `ind`
            status_ind
                Simulation end status of repetition no. `ind`
        """
        if ind in self:
            x_ind = self.x_list[ind]
            t_ind = self.t_list[ind]
            status_ind = self.status_list[ind]
        else:
            raise IndexError(f"{ind} out of bounds")
        return x_ind, t_ind, status_ind

    @property
    def final(self) -> Tuple[np.ndarray, np.ndarray]:
        """
            Returns the final times and states of the system in the simulations

            Returns
            -------
            Tuple[np.ndarray, np.ndarray]
                The final times and states of the sytem
        """
        final_times = np.array([v[1][-1] for v in self])
        final_states = np.array([v[0][-1, :] for v in self])
        return final_times, final_states

    def get_state(self, t: float) -> List[np.ndarray]:
        """
            Returns the states of the system at time point t.

            Parameters
            ----------
            t: float
                Time point at which states are wanted.

            Returns
            -------
            List[np.ndarray]
                The states of the system at `t` for all repetitions.

            Raises
            ------
            UserWarning
                If simulation ends before `t` but system does not reach
                extinction.
        """
        states: List[np.ndarray] = []
        for x_array, t_array, s in self:
            ind = np.searchsorted(t_array, t)
            ind = ind - 1 if ind > 0 else ind
            if ind == len(t_array) - 1:
                states.append(x_array[-1, :])
                if s != 3:
                    warn(f"Simulation ended before {t}, returning last state.")
            else:
                x_interp = np.zeros(x_array.shape[1])
                if self.algorithm != "direct":
                    for ind2 in range(x_array.shape[1]):
                        x_interp[ind2] = np.interp(
                            t,
                            [t_array[ind], t_array[ind + 1]],
                            [x_array[ind, ind2], x_array[ind + 1, ind2]],
                        )
                    states.append(x_interp)
                else:
                    states.append(x_array[ind, :])
        return states
