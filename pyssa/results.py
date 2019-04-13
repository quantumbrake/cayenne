"""
    Module that defines the `Results` class
"""

from collections.abc import Collection
from typing import List, Tuple

import numpy as np


class Results(Collection):
    """
        A class that stores simulation results and provides methods to access them

        Parameters
        ----------
        t_list : List[float]
        x_list : List[np.ndarray]
        status_list : List[int]
        algorithm : str
        seed: List[int]

        Other Parameters
        ----------------

        Attributes
        ----------
    """

    def __init__(
        self,
        t_list: List[np.ndarray],
        x_list: List[np.ndarray],
        status_list: List[int],
        algorithm: str,
        seed: List[int],
        **kwargs,
    ) -> None:
        self.x_list = x_list
        self.t_list = t_list
        self.status_list = status_list
        self.algorithm = algorithm
        self.seed = seed
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
            == len(self.seed)
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
        return (
            f"<Results n_rep={len(self)} algorithm={self.algorithm} seed={self.seed}>"
        )

    def __str__(self) -> str:
        return self.__repr__()

    def __iter__(self):
        return zip(self.x_list, self.t_list, self.status_list)

    def __len__(self):
        return len(self.x_list)

    def __contains__(self, ind):
        if ind < len(self):
            return True
        else:
            return False

    def __getitem__(self, ind: int):
        if ind in self:
            return self.x_list[ind], self.t_list[ind], self.status_list[ind]
        else:
            raise IndexError(f"{ind} out of bounds")

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
            t : float
                Time point at which states are wanted.

            Returns
            -------
            List[np.ndarray]
                The states of the system at `t` for all repetitions.
        """
        states: List[np.ndarray] = []
        for x_array, t_array, _ in self:
            ind = np.searchsorted(t_array, t)
            ind = ind - 1 if ind > 0 else ind
            x_interp = np.zeros(x_array.shape[1])
            for ind2 in range(x_array.shape[1]):
                x_interp[ind2] = np.interp(
                    t,
                    [t_array[ind], t_array[ind + 1]],
                    [x_array[ind, ind2], x_array[ind + 1, ind2]],
                )
            states.append(x_interp)
        return states
